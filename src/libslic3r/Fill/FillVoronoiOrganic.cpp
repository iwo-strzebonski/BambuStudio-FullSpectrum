// Voronoi Organic infill.
//
// Generates a Voronoi tessellation whose seed-point density is driven by the
// distance to the nearest boundary wall.  Seeds cluster near walls (high
// simulated bending-stress zone) and spread apart in the interior (low-stress
// bulk).  The resulting irregular cell network is clipped to the ExPolygon and
// output as printable line segments.
//
// Coordinate conventions:
//   * All geometry for the Voronoi diagram is kept in UNSCALED millimetres
//     (double) to avoid precision problems in boost::polygon.
//   * We only re-scale to Slic3r scaled coords when constructing Points for
//     the final Polylines.

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include "FillVoronoiOrganic.hpp"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146)  // unary minus on unsigned type (boost)
#endif
#include <boost/polygon/voronoi.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <cmath>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>

namespace Slic3r {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Squared distance from point P to line-segment AB (all in unscaled mm).
static double point_segment_dist_sq(const Vec2d &p, const Vec2d &a, const Vec2d &b)
{
    Vec2d ab = b - a;
    double len_sq = ab.squaredNorm();
    if (len_sq < 1e-14)
        return (p - a).squaredNorm();
    double t = std::clamp((p - a).dot(ab) / len_sq, 0.0, 1.0);
    Vec2d proj = a + ab * t;
    return (p - proj).squaredNorm();
}

// Minimum squared distance from P to any edge of any ring in the ExPolygon.
// Input/output in unscaled mm (because ExPolygon coords are scaled, we
// unscale internally).
static double point_to_expolygon_dist_sq(const Vec2d &p_mm, const ExPolygon &exp)
{
    double dmin = std::numeric_limits<double>::max();

    auto process_ring = [&](const Polygon &ring) {
        const Points &pts = ring.points;
        size_t n = pts.size();
        for (size_t i = 0; i < n; ++i) {
            Vec2d a = unscale(pts[i]);
            Vec2d b = unscale(pts[(i + 1) % n]);
            double d = point_segment_dist_sq(p_mm, a, b);
            if (d < dmin)
                dmin = d;
        }
    };

    process_ring(exp.contour);
    for (const Polygon &hole : exp.holes)
        process_ring(hole);

    return dmin;
}

// ---------------------------------------------------------------------------
// Spatial-hash cell for Poisson-disk rejection sampling.
// Key encodes the 2-D integer cell index into a single 64-bit value.
// ---------------------------------------------------------------------------
static inline int64_t spatial_key(int ix, int iy)
{
    // clamp to 20-bit range so the pack fits in 40 bits and is collision-free
    return (int64_t(ix & 0xFFFFF) << 20) | (int64_t(iy & 0xFFFFF));
}

// ---------------------------------------------------------------------------
// Seed generation
// ---------------------------------------------------------------------------

// Returns a vector of Voronoi seed points (unscaled mm) distributed so that
// their local spacing reflects the simulated stress field:
//   - dense near the boundary walls (distance ≈ 0)
//   - sparse in the interior (distance ≫ falloff)
// The z argument introduces per-layer jitter via low-frequency trigonometry.
static std::vector<Vec2d> generate_stress_seeds(
    const ExPolygon &expolygon,
    double           base_spacing,   // unscaled mm
    double           z)              // layer height (mm), for jitter
{
    const double dense_sp  = base_spacing * FillVoronoiOrganic::DenseFactor;
    const double sparse_sp = base_spacing * FillVoronoiOrganic::SparseFactor;
    const double falloff   = base_spacing * FillVoronoiOrganic::FalloffFactor;

    // Per-layer jitter offsets — vary slowly with Z for an organic look.
    const double jitter_x = std::cos(z * 3.7) * dense_sp * FillVoronoiOrganic::JitterAmplitude;
    const double jitter_y = std::sin(z * 2.9) * dense_sp * FillVoronoiOrganic::JitterAmplitude;

    BoundingBoxf bb = unscaled(get_extents(expolygon));

    // Grid step = dense_sp (worst-case) so we don't skip any candidate.
    const double step = dense_sp;
    if (step < 1e-6)
        return {};

    // spatial hash: cell_size = dense_sp / sqrt(2) ensures only same/adjacent
    // cells need to be checked for the acceptance radius.
    const double cell_size = dense_sp * 0.7072;

    // Store accepted seed coords (unscaled mm).
    // hash: cell_key → list of seeds in that cell
    std::unordered_map<int64_t, std::vector<Vec2d>> spatial;
    std::vector<Vec2d> seeds;

    // Helper: query whether a candidate is too close to any existing seed
    // within required_spacing.
    auto has_neighbor = [&](const Vec2d &cand, double req_sp) -> bool {
        int cx = (int)std::floor((cand.x() - bb.min.x()) / cell_size);
        int cy = (int)std::floor((cand.y() - bb.min.y()) / cell_size);
        int r  = (int)std::ceil(req_sp / cell_size) + 1;
        for (int dx = -r; dx <= r; ++dx) {
            for (int dy = -r; dy <= r; ++dy) {
                auto it = spatial.find(spatial_key(cx + dx, cy + dy));
                if (it == spatial.end())
                    continue;
                for (const Vec2d &s : it->second)
                    if ((cand - s).squaredNorm() < sqr(req_sp))
                        return true;
            }
        }
        return false;
    };

    auto insert_seed = [&](const Vec2d &pt) {
        int cx = (int)std::floor((pt.x() - bb.min.x()) / cell_size);
        int cy = (int)std::floor((pt.y() - bb.min.y()) / cell_size);
        spatial[spatial_key(cx, cy)].push_back(pt);
        seeds.push_back(pt);
    };

    // Iterate over candidate grid positions.
    for (double gx = bb.min.x(); gx <= bb.max.x() + step; gx += step) {
        for (double gy = bb.min.y(); gy <= bb.max.y() + step; gy += step) {
            // Apply Z-based jitter so adjacent layers look different.
            Vec2d cand(gx + jitter_x, gy + jitter_y);

            // Only accept points inside the shape.
            Point scaled_cand(scale_(cand.x()), scale_(cand.y()));
            if (!expolygon.contains(scaled_cand))
                continue;

            // Compute distance-derived spacing at this point.
            double dist_sq = point_to_expolygon_dist_sq(cand, expolygon);
            double dist    = std::sqrt(dist_sq);
            // t=0 → wall, t=1 → interior
            double t       = std::min(dist / falloff, 1.0);
            t              = t * t;  // quadratic falloff for crisper edge zone
            double req_sp  = dense_sp + t * (sparse_sp - dense_sp);

            if (!has_neighbor(cand, req_sp))
                insert_seed(cand);
        }
    }

    return seeds;
}

// ---------------------------------------------------------------------------
// Infinite Voronoi edge clipping
// ---------------------------------------------------------------------------

// Given the two sites defining an infinite Voronoi edge, and one finite
// vertex on it (if any), produce a clipped Vec2d endpoint within the bbox.
// If the edge has no finite vertex at all (both infinite), the midpoint of
// the two sites is the foot of the perpendicular.
static Vec2d clip_ray_to_bbox(Vec2d foot, Vec2d dir, const BoundingBoxf &bbox, bool negative)
{
    if (negative)
        dir = -dir;

    // Normalise direction
    double len = dir.norm();
    if (len < 1e-14)
        return foot;
    dir /= len;

    // Pad the bbox slightly so clipped edges definitely exit the shape.
    double large = std::max(bbox.max.x() - bbox.min.x(),
                            bbox.max.y() - bbox.min.y()) * 2.0 + 10.0;

    // Walk large distance along direction; the actual clip will come from
    // intersection_pl with the ExPolygon, so we just need a far-enough point.
    return foot + dir * large;
}

// ---------------------------------------------------------------------------
// Main fill method
// ---------------------------------------------------------------------------

void FillVoronoiOrganic::_fill_surface_single(
    const FillParams              &params,
    unsigned int                   /*thickness_layers*/,
    const std::pair<float, Point> & /*direction*/,
    ExPolygon                       expolygon,
    Polylines                      &polylines_out)
{
    // base_spacing in unscaled mm — respects infill density.
    const double base_spacing_mm = this->spacing / params.density;

    // Generate stress-following seeds (unscaled mm coords).
    std::vector<Vec2d> seeds = generate_stress_seeds(expolygon, base_spacing_mm, this->z);

    // Need at least 3 seeds for a meaningful Voronoi diagram.
    if (seeds.size() < 3)
        return;

    // Build boost::polygon Voronoi diagram from point sites.
    using BoostPoint = boost::polygon::point_data<double>;
    std::vector<BoostPoint> bpoints;
    bpoints.reserve(seeds.size());
    for (const Vec2d &s : seeds)
        bpoints.emplace_back(s.x(), s.y());

    boost::polygon::voronoi_diagram<double> vd;
    boost::polygon::construct_voronoi(bpoints.begin(), bpoints.end(), &vd);

    // BoundingBoxf for infinite-edge clipping (inflated slightly).
    BoundingBoxf bb = unscaled(get_extents(expolygon));
    {
        double pad = base_spacing_mm * 2.0;
        bb.min -= Vec2d(pad, pad);
        bb.max += Vec2d(pad, pad);
    }

    // Collect Voronoi edge segments as Polylines (unscaled mm → scaled Points).
    // Each undirected edge appears as two half-edges; we only process the one
    // where cell source_index < twin source_index to avoid duplicates.
    Polylines raw_edges;
    raw_edges.reserve(vd.num_edges() / 2);

    for (const auto &edge : vd.edges()) {
        // Skip secondary (degenerate/infinite internal) edges.
        if (!edge.is_primary())
            continue;

        // Deduplicate twin pairs.
        const auto *twin = edge.twin();
        if (edge.cell()->source_index() >= twin->cell()->source_index())
            continue;

        // Site midpoint and perpendicular — needed for infinite edge direction.
        const Vec2d &s0 = seeds[edge.cell()->source_index()];
        const Vec2d &s1 = seeds[twin->cell()->source_index()];
        Vec2d        diff(s1.x() - s0.x(), s1.y() - s0.y());
        Vec2d        perp(-diff.y(), diff.x());  // perpendicular to s0→s1
        Vec2d        foot((s0.x() + s1.x()) * 0.5, (s0.y() + s1.y()) * 0.5);

        Vec2d p0, p1;

        const auto *v0 = edge.vertex0();
        const auto *v1 = edge.vertex1();

        if (v0 && v1) {
            // Fully finite edge.
            p0 = Vec2d(v0->x(), v0->y());
            p1 = Vec2d(v1->x(), v1->y());
        } else if (v0) {
            // v1 is at infinity — clip from v0 outward.
            p0 = Vec2d(v0->x(), v0->y());
            p1 = clip_ray_to_bbox(p0, perp, bb, false);
        } else if (v1) {
            // v0 is at infinity — clip from v1 inward.
            p1 = Vec2d(v1->x(), v1->y());
            p0 = clip_ray_to_bbox(p1, perp, bb, true);
        } else {
            // Both vertices at infinity — extend from midpoint in both dirs.
            p0 = clip_ray_to_bbox(foot, perp, bb, true);
            p1 = clip_ray_to_bbox(foot, perp, bb, false);
        }

        // Convert to scaled Slic3r coords and store.
        Polyline pl;
        pl.points.emplace_back(scale_(p0.x()), scale_(p0.y()));
        pl.points.emplace_back(scale_(p1.x()), scale_(p1.y()));
        raw_edges.push_back(std::move(pl));
    }

    if (raw_edges.empty())
        return;

    // Clip all edges to the ExPolygon boundary.
    Polylines clipped = intersection_pl(raw_edges, expolygon);

    // Discard degenerate segments shorter than half the minimum cell size.
    const coord_t min_len = scale_(base_spacing_mm * FillVoronoiOrganic::DenseFactor * 0.5);
    clipped.erase(
        std::remove_if(clipped.begin(), clipped.end(),
            [min_len](const Polyline &pl) { return pl.length() < min_len; }),
        clipped.end());

    if (clipped.empty())
        return;

    // Sort segments by travel proximity then emit.
    append(polylines_out, chain_polylines(std::move(clipped)));
}

} // namespace Slic3r
