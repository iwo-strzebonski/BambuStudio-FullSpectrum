// Cellular Foam infill.
//
// Generates a Voronoi tessellation whose seed points are distributed
// stochastically at two scales to produce irregular polygon cells of varying
// sizes — mimicking a cross-section through closed-cell foam.
//
// Coordinate conventions:
//   * Seed placement and Voronoi construction use UNSCALED millimetres (double).
//   * We only convert to Slic3r scaled int coords when building the output
//     Polylines, matching the convention used by FillVoronoiOrganic.

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include "FillCellularFoam.hpp"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146)  // unary minus on unsigned type (boost)
#endif
#include <boost/polygon/voronoi.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>

namespace Slic3r {

// ---------------------------------------------------------------------------
// PCG32 — minimal Permuted Congruential Generator
// Fast, high-quality 32-bit PRNG; seeded per-layer for reproducible patterns.
// Reference: O'Neill (2014) https://www.pcg-random.org/
// ---------------------------------------------------------------------------
struct Pcg32 {
    uint64_t state;
    uint64_t inc;

    explicit Pcg32(uint64_t seed, uint64_t seq = 1442695040888963407ULL)
    {
        state = 0;
        inc   = (seq << 1u) | 1u;
        step();
        state += seed;
        step();
    }

    uint32_t step()
    {
        uint64_t old = state;
        state = old * 6364136223846793005ULL + inc;
        uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
        uint32_t rot        = (uint32_t)(old >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-(int32_t)rot) & 31));
    }

    // Returns a float in [0, 1).
    float nextf() { return (float)(step() >> 8) * (1.0f / (float)(1 << 24)); }
};

// ---------------------------------------------------------------------------
// Spatial-hash utilities — identical layout to FillVoronoiOrganic helpers
// to share the same proven spatial-rejection logic.
// ---------------------------------------------------------------------------
static inline int64_t foam_spatial_key(int ix, int iy)
{
    return (int64_t(ix & 0xFFFFF) << 20) | int64_t(iy & 0xFFFFF);
}

// ---------------------------------------------------------------------------
// Multi-scale stochastic seed generation
// ---------------------------------------------------------------------------
// Returns seed points (unscaled mm) representing foam cell centres sampled
// at two spatial scales.  The coarse pass produces large cells; the fine pass
// fills the interstices with smaller cells.
//
// The PCG32 is seeded from layer_id and the bounding-box origin so every layer
// and every print object gets a unique yet reproducible random pattern.
static std::vector<Vec2d> generate_foam_seeds(
    const ExPolygon &expolygon,
    double           base_spacing,  // unscaled mm, controls overall cell density
    unsigned int     layer_id)
{
    BoundingBoxf bb = unscaled(get_extents(expolygon));

    // --- PRNG seeded from layer index and bounding-box position ---------------
    // Multiply by large primes so small coordinate values still spread well
    // across the seed space.
    uint64_t prng_seed =
        (uint64_t)layer_id * 2654435761ULL
        ^ (uint64_t)(std::abs(bb.min.x()) * 997.0)  * 3141592653ULL
        ^ (uint64_t)(std::abs(bb.min.y()) * 997.0)  * 2718281829ULL;
    Pcg32 rng(prng_seed);

    // --- Pass parameters -----------------------------------------------------
    const double coarse_step    = base_spacing * FillCellularFoam::CoarseScale;
    const double fine_step      = base_spacing * FillCellularFoam::FineScale;
    const double min_sep        = base_spacing * FillCellularFoam::MinSepFrac;

    // Spatial-hash grid cell size: min_sep / sqrt(2) ensures that only
    // same-cell and 8-adjacent cells need to be checked.
    const double cell_size = min_sep * 0.7072;
    if (cell_size < 1e-6)
        return {};

    std::unordered_map<int64_t, std::vector<Vec2d>> spatial;
    std::vector<Vec2d> seeds;

    // --- Spatial-hash helpers ------------------------------------------------
    auto has_neighbor = [&](const Vec2d &cand) -> bool {
        int cx = (int)std::floor((cand.x() - bb.min.x()) / cell_size);
        int cy = (int)std::floor((cand.y() - bb.min.y()) / cell_size);
        int r  = (int)std::ceil(min_sep / cell_size) + 1;
        for (int dx = -r; dx <= r; ++dx) {
            for (int dy = -r; dy <= r; ++dy) {
                auto it = spatial.find(foam_spatial_key(cx + dx, cy + dy));
                if (it == spatial.end())
                    continue;
                for (const Vec2d &s : it->second)
                    if ((cand - s).squaredNorm() < min_sep * min_sep)
                        return true;
            }
        }
        return false;
    };

    auto insert_seed = [&](const Vec2d &pt) {
        int cx = (int)std::floor((pt.x() - bb.min.x()) / cell_size);
        int cy = (int)std::floor((pt.y() - bb.min.y()) / cell_size);
        spatial[foam_spatial_key(cx, cy)].push_back(pt);
        seeds.push_back(pt);
    };

    // =========================================================================
    // Pass 1 — Coarse seeds (large bubbles)
    // =========================================================================
    // Step slightly outside the bounding box so cells that straddle the edges
    // still get clipped polygon walls rather than abruptly terminated edges.
    const double coarse_accept = (float)FillCellularFoam::CoarseAcceptProb;
    const double coarse_jitter = coarse_step * FillCellularFoam::JitterFrac;

    for (double gx = bb.min.x() - coarse_step; gx <= bb.max.x() + coarse_step; gx += coarse_step) {
        for (double gy = bb.min.y() - coarse_step; gy <= bb.max.y() + coarse_step; gy += coarse_step) {
            // Stochastic acceptance — creates the size variation in large cells.
            if (rng.nextf() > coarse_accept)
                continue;

            // Random jitter: ± coarse_jitter in each axis.
            Vec2d cand(
                gx + (rng.nextf() - 0.5f) * 2.0 * coarse_jitter,
                gy + (rng.nextf() - 0.5f) * 2.0 * coarse_jitter);

            // Accept only seeds inside the boundary.
            Point sc(scale_(cand.x()), scale_(cand.y()));
            if (!expolygon.contains(sc))
                continue;

            if (!has_neighbor(cand))
                insert_seed(cand);
        }
    }

    // =========================================================================
    // Pass 2 — Fine seeds (small bubbles that fill the interstices)
    // =========================================================================
    const double fine_accept = (float)FillCellularFoam::FineAcceptProb;
    const double fine_jitter = fine_step * FillCellularFoam::JitterFrac;

    for (double gx = bb.min.x() - fine_step; gx <= bb.max.x() + fine_step; gx += fine_step) {
        for (double gy = bb.min.y() - fine_step; gy <= bb.max.y() + fine_step; gy += fine_step) {
            if (rng.nextf() > fine_accept)
                continue;

            Vec2d cand(
                gx + (rng.nextf() - 0.5f) * 2.0 * fine_jitter,
                gy + (rng.nextf() - 0.5f) * 2.0 * fine_jitter);

            Point sc(scale_(cand.x()), scale_(cand.y()));
            if (!expolygon.contains(sc))
                continue;

            if (!has_neighbor(cand))
                insert_seed(cand);
        }
    }

    return seeds;
}

// ---------------------------------------------------------------------------
// Clip an infinite Voronoi ray to a bounding box, identical to the helper
// in FillVoronoiOrganic (duplicated to keep files self-contained).
// ---------------------------------------------------------------------------
static Vec2d foam_clip_ray(Vec2d foot, Vec2d dir, const BoundingBoxf &bbox, bool negative)
{
    if (negative)
        dir = -dir;
    double len = dir.norm();
    if (len < 1e-14)
        return foot;
    dir /= len;
    double large = (std::max(bbox.max.x() - bbox.min.x(),
                             bbox.max.y() - bbox.min.y()) * 2.0) + 10.0;
    return foot + dir * large;
}

// ---------------------------------------------------------------------------
// Main fill entry point
// ---------------------------------------------------------------------------
void FillCellularFoam::_fill_surface_single(
    const FillParams              &params,
    unsigned int                   /*thickness_layers*/,
    const std::pair<float, Point> & /*direction*/,
    ExPolygon                       expolygon,
    Polylines                      &polylines_out)
{
    // Base cell spacing in unscaled mm — larger spacing == fewer cells == lower density.
    const double base_spacing_mm = this->spacing / params.density;

    // Generate stochastic multi-scale foam seeds.
    std::vector<Vec2d> seeds = generate_foam_seeds(expolygon, base_spacing_mm, this->layer_id);

    // Voronoi requires at least 3 sites.
    if (seeds.size() < 3)
        return;

    // --- Build Voronoi diagram -----------------------------------------------
    using BoostPoint = boost::polygon::point_data<double>;
    std::vector<BoostPoint> bpoints;
    bpoints.reserve(seeds.size());
    for (const Vec2d &s : seeds)
        bpoints.emplace_back(s.x(), s.y());

    boost::polygon::voronoi_diagram<double> vd;
    boost::polygon::construct_voronoi(bpoints.begin(), bpoints.end(), &vd);

    // Inflated bounding box for clipping infinite Voronoi rays.
    BoundingBoxf bb = unscaled(get_extents(expolygon));
    {
        double pad = base_spacing_mm * 2.0;
        bb.min -= Vec2d(pad, pad);
        bb.max += Vec2d(pad, pad);
    }

    // --- Extract & convert edges ---------------------------------------------
    // Process each undirected edge exactly once: take the half-edge whose
    // cell source_index is strictly less than its twin's source_index.
    Polylines raw_edges;
    raw_edges.reserve(vd.num_edges() / 2);

    for (const auto &edge : vd.edges()) {
        if (!edge.is_primary())
            continue;

        const auto *twin = edge.twin();
        if (edge.cell()->source_index() >= twin->cell()->source_index())
            continue;

        const Vec2d &s0 = seeds[edge.cell()->source_index()];
        const Vec2d &s1 = seeds[twin->cell()->source_index()];

        // Perpendicular bisector direction (needed for infinite edges).
        Vec2d diff(s1.x() - s0.x(), s1.y() - s0.y());
        Vec2d perp(-diff.y(), diff.x());
        Vec2d foot((s0.x() + s1.x()) * 0.5, (s0.y() + s1.y()) * 0.5);

        Vec2d p0, p1;
        const auto *v0 = edge.vertex0();
        const auto *v1 = edge.vertex1();

        if (v0 && v1) {
            p0 = Vec2d(v0->x(), v0->y());
            p1 = Vec2d(v1->x(), v1->y());
        } else if (v0) {
            p0 = Vec2d(v0->x(), v0->y());
            p1 = foam_clip_ray(p0, perp, bb, false);
        } else if (v1) {
            p1 = Vec2d(v1->x(), v1->y());
            p0 = foam_clip_ray(p1, perp, bb, true);
        } else {
            p0 = foam_clip_ray(foot, perp, bb, true);
            p1 = foam_clip_ray(foot, perp, bb, false);
        }

        Polyline pl;
        pl.points.emplace_back(scale_(p0.x()), scale_(p0.y()));
        pl.points.emplace_back(scale_(p1.x()), scale_(p1.y()));
        raw_edges.push_back(std::move(pl));
    }

    if (raw_edges.empty())
        return;

    // --- Clip to ExPolygon ---------------------------------------------------
    Polylines clipped = intersection_pl(raw_edges, expolygon);

    // Discard segments that are too short to print meaningfully.
    const coord_t min_len = scale_(base_spacing_mm * FillCellularFoam::MinSegFrac);
    clipped.erase(
        std::remove_if(clipped.begin(), clipped.end(),
            [min_len](const Polyline &pl) { return pl.length() < min_len; }),
        clipped.end());

    if (clipped.empty())
        return;

    // --- Sort for travel efficiency and emit ---------------------------------
    append(polylines_out, chain_polylines(std::move(clipped)));
}

} // namespace Slic3r
