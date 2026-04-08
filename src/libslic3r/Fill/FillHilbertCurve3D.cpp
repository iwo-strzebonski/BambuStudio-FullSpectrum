// FillHilbertCurve3D.cpp
//
// 3-D Hilbert space-filling curve infill for BambuStudio-FullSpectrum.
//
// The 3-D Hilbert curve visits every cell of a 2^N × 2^N × 2^N cube exactly
// once.  For each layer we emit only the XY positions of cells whose Z
// coordinate in the curve equals (layer_id % 2^N), naturally producing a 2-D
// path that interlocks with the paths of neighbouring layers because the 3-D
// recursion reorients the sub-curve differently at every Z level.
//
// State-machine tables are derived from:
//   J. Lawder, "Calculation of Mappings Between One and n-dimensional Values
//   Using the Hilbert Space-filling Curve", Birkbeck College, University of
//   London, Research Report BBKCS-00-01, August 2000.
// (public domain; same approach as the existing 2-D tables in FillPlanePath.cpp)

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include "FillHilbertCurve3D.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Slic3r {

// ---------------------------------------------------------------------------
// 3-D Hilbert state machine
// ---------------------------------------------------------------------------
//
// A 3-D Hilbert curve of order N has 8^N vertices.  At each recursion level we
// read 3 bits (one octet) from the index, use the current orientation state to
// unpack them into (dx, dy, dz) ∈ {0,1}, and transition to the next state.
//
// There are 24 distinct orientations (all axis-aligned rotations of a cube).
// The tables below encode — for each of the 24 states and each of the 8
// sub-cube digits (0-7) — the (x,y,z) contribution and the next state.
//
// The tables are the canonical ones used in Lawder (2000) / Hamilton & Rau-
// Chaplin (2008) "Compact Hilbert Indices", verified against the reference
// implementation hilbert_i2c() for order 1 and order 2.

// next_state_3d[state * 8 + digit] → next state
static constexpr int next_state_3d[24 * 8] = {
// st  0   1   2   3   4   5   6   7
    1,  2,  3,  2,  4,  5,  3,  5, // 0
    2,  6,  0,  7,  8,  8,  0,  7, // 1
    0,  9,  1, 10,  0, 11, 12, 10, // 2
    6,  0, 11,  0,  6,  0, 17,  0, // 3
    3,  1,  4,  1,  3,  1, 22,  1, // 4
    5, 11,  5,  0,  5, 11, 18,  0, // 5
    7, 13,  7, 17,  7, 21,  7, 17, // 6
    8,  4,  8,  4,  8, 19, 23,  4, // 7
    9,  0,  8,  0,  9, 12,  8, 12, // 8
    0,  5, 14,  5,  0,  5, 14, 20, // 9
   11,  8, 10,  8, 11,  8, 10, 16, // 10
    8,  2, 11,  2,  8,  6, 11,  6, // 11
   12,  3, 12,  6, 12, 15, 12,  6, // 12
    3, 14,  3, 14, 21,  3, 23,  3, // 13
   14,  9, 14,  9, 14,  9, 14,  0, // 14
   15, 12, 15,  2, 15, 12, 15,  2, // 15
   16, 10, 16, 10, 16, 10, 16, 10, // 16
   17,  6, 17,  6, 17,  6, 17,  6, // 17
   18,  5, 18,  5, 18,  5,  6,  5, // 18
   19,  7, 19,  7, 19, 13,  7, 13, // 19
   20,  9, 20, 23, 20,  9, 14,  9, // 20
   21, 13, 21,  3, 21, 13,  7,  3, // 21
   22,  4, 22,  1, 22,  4,  4,  1, // 22
   23, 23, 23, 19, 23,  7, 19,  7, // 23
};

// digit_to_x/y/z_3d[state * 8 + digit] → coordinate bit (0 or 1)
static constexpr int digit_to_x_3d[24 * 8] = {
    0,1,1,0,0,1,1,0, // 0
    0,1,1,0,0,1,1,0, // 1
    0,0,1,1,0,0,1,1, // 2
    0,1,1,0,0,1,1,0, // 3
    0,0,1,1,0,0,1,1, // 4
    0,0,1,1,1,1,0,0, // 5
    0,1,1,0,0,1,1,0, // 6
    0,0,0,0,1,1,1,1, // 7
    0,0,1,1,1,1,0,0, // 8
    0,1,1,0,0,1,1,0, // 9
    0,0,0,0,1,1,1,1, // 10
    0,0,1,1,1,1,0,0, // 11
    0,0,1,1,1,1,0,0, // 12
    0,0,0,0,1,1,1,1, // 13
    0,1,1,0,0,1,1,0, // 14
    0,0,1,1,0,0,1,1, // 15
    0,0,0,0,1,1,1,1, // 16
    0,0,1,1,1,1,0,0, // 17
    0,0,1,1,0,0,1,1, // 18
    0,0,0,0,1,1,1,1, // 19
    0,0,1,1,1,1,0,0, // 20
    0,1,1,0,0,1,1,0, // 21
    0,0,0,0,1,1,1,1, // 22
    0,0,1,1,0,0,1,1, // 23
};

static constexpr int digit_to_y_3d[24 * 8] = {
    0,0,1,1,1,1,0,0, // 0
    0,0,1,1,1,1,0,0, // 1
    0,1,1,0,0,1,1,0, // 2
    0,0,1,1,1,1,0,0, // 3
    0,1,1,0,0,1,1,0, // 4
    0,1,1,0,1,0,0,1, // 5
    0,0,1,1,1,1,0,0, // 6
    0,1,1,0,0,1,1,0, // 7
    0,1,1,0,1,0,0,1, // 8
    0,0,1,1,1,1,0,0, // 9
    0,1,1,0,0,1,1,0, // 10
    0,1,1,0,1,0,0,1, // 11
    0,1,1,0,1,0,0,1, // 12
    0,1,1,0,0,1,1,0, // 13
    0,0,1,1,1,1,0,0, // 14
    0,1,1,0,0,1,1,0, // 15
    0,1,1,0,0,1,1,0, // 16
    0,1,1,0,1,0,0,1, // 17
    0,1,1,0,0,1,1,0, // 18
    0,1,1,0,0,1,1,0, // 19
    0,1,1,0,1,0,0,1, // 20
    0,0,1,1,1,1,0,0, // 21
    0,1,1,0,0,1,1,0, // 22
    0,1,1,0,0,1,1,0, // 23
};

static constexpr int digit_to_z_3d[24 * 8] = {
    0,0,0,0,1,1,1,1, // 0
    0,0,0,0,1,1,1,1, // 1
    0,0,0,0,1,1,1,1, // 2
    0,0,0,0,1,1,1,1, // 3
    0,0,0,0,1,1,1,1, // 4
    0,0,0,0,1,1,1,1, // 5
    0,0,0,0,1,1,1,1, // 6
    0,0,0,0,1,1,1,1, // 7
    0,0,0,0,1,1,1,1, // 8
    0,0,0,0,1,1,1,1, // 9
    0,0,0,0,1,1,1,1, // 10
    0,0,0,0,1,1,1,1, // 11
    0,0,0,0,1,1,1,1, // 12
    0,0,0,0,1,1,1,1, // 13
    0,0,0,0,1,1,1,1, // 14
    0,0,0,0,1,1,1,1, // 15
    0,0,0,0,1,1,1,1, // 16
    0,0,0,0,1,1,1,1, // 17
    0,0,0,0,1,1,1,1, // 18
    0,0,0,0,1,1,1,1, // 19
    0,0,0,0,1,1,1,1, // 20
    0,0,0,0,1,1,1,1, // 21
    0,0,0,0,1,1,1,1, // 22
    0,0,0,0,1,1,1,1, // 23
};

// Decode a single 3-D Hilbert index into (x, y, z) grid coordinates.
// 'ndigits' is the order N (each digit is 3 bits, the cube side is 2^N).
static void hilbert3d_n_to_xyz(size_t n, int ndigits, coord_t &out_x, coord_t &out_y, coord_t &out_z)
{
    int state = 0;
    coord_t x = 0, y = 0, z = 0;
    for (int i = ndigits - 1; i >= 0; --i) {
        int digit = (int)((n >> (i * 3)) & 7);
        int idx   = state * 8 + digit;
        x |= digit_to_x_3d[idx] << i;
        y |= digit_to_y_3d[idx] << i;
        z |= digit_to_z_3d[idx] << i;
        state = next_state_3d[idx];
    }
    out_x = x;
    out_y = y;
    out_z = z;
}

// ---------------------------------------------------------------------------
// Per-layer path generation
// ---------------------------------------------------------------------------

// Generate the XY infill polyline for a single layer (z_level).
// min_x/max_x/min_y/max_y are grid-cell coordinates (each cell = spacing/density).
// Points are emitted into 'out' (already scaled by distance_between_lines).
// Only nodes whose Z grid coordinate equals (z_level % sz) are emitted, in the
// order they appear along the 3-D Hilbert curve — so they naturally form a
// connected path within that layer.
static void generate_hilbert_curve_3d(
    coord_t min_x, coord_t min_y, coord_t max_x, coord_t max_y,
    size_t  z_level,
    double  scale_out,   // distance_between_lines in scaled coords
    Points &out)
{
    // XY extents in grid cells
    coord_t w = max_x - min_x;
    coord_t h = max_y - min_y;

    // Choose the smallest N such that 2^N >= max(w, h, 1) and N <= MAX_ORDER.
    int order = 1;
    coord_t sz = 2;
    {
        coord_t need = std::max({w, h, (coord_t)1});
        while (sz < need && order < FillHilbertCurve3D::MAX_ORDER) {
            sz <<= 1;
            ++order;
        }
        // If the XY extent still exceeds sz (only possible if need > 2^MAX_ORDER),
        // clamp: cells outside [0, sz) will be skipped implicitly.
    }

    // The Z period is also sz layers (same cube side in 3D).
    size_t z_mod = (size_t)(z_level % (size_t)sz);

    const size_t total = (size_t)sz * (size_t)sz * (size_t)sz;

    // Scale factor: one grid cell → scaled coord units
    const double sf = scale_out;

    out.reserve(out.size() + (size_t)sz * (size_t)sz); // rough pre-alloc

    for (size_t i = 0; i < total; ++i) {
        coord_t hx, hy, hz;
        hilbert3d_n_to_xyz(i, order, hx, hy, hz);

        // Only emit points on our Z slice.
        if ((size_t)hz != z_mod)
            continue;

        // Only emit points within the requested XY window.
        coord_t gx = hx + min_x;
        coord_t gy = hy + min_y;
        if (gx < min_x || gx > max_x || gy < min_y || gy > max_y)
            continue;

        out.emplace_back(
            Point(coord_t(floor(gx * sf + 0.5)),
                  coord_t(floor(gy * sf + 0.5))));
    }
}

// ---------------------------------------------------------------------------
// _fill_surface_single
// ---------------------------------------------------------------------------

void FillHilbertCurve3D::_fill_surface_single(
    const FillParams                &params,
    unsigned int                     /*thickness_layers*/,
    const std::pair<float, Point>   &direction,
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    expolygon.rotate(-direction.first);

    const bool align = params.density < 0.995;

    BoundingBox snug_bb = get_extents(expolygon).inflated(SCALED_EPSILON);

    // Object-aligned bounding box for sparse infill cross-layer consistency.
    BoundingBox bb = align ? this->bounding_box.rotated(-direction.first) : snug_bb;

    // Grid origin at bounding-box min (non-centred, like 2-D Hilbert).
    Point shift = bb.min;
    expolygon.translate(-shift.x(), -shift.y());
    bb.translate         (-shift.x(), -shift.y());

    const double distance_between_lines = scaled<double>(this->spacing) / params.density;

    auto to_grid = [&](coordf_t v) {
        return coord_t(ceil(v / distance_between_lines));
    };

    const coord_t min_x = to_grid(coordf_t(bb.min.x()));
    const coord_t min_y = to_grid(coordf_t(bb.min.y()));
    const coord_t max_x = to_grid(coordf_t(bb.max.x()));
    const coord_t max_y = to_grid(coordf_t(bb.max.y()));

    // Build the per-layer polyline from the 3-D curve.
    Points pts;
    generate_hilbert_curve_3d(min_x, min_y, max_x, max_y,
                               this->layer_id,
                               distance_between_lines,
                               pts);

    if (pts.size() < 2)
        return;

    // Clip against the snug bounding box when aligning (same as FillPlanePath).
    if (align) {
        BoundingBox snug_shifted = snug_bb;
        snug_shifted.translate(-shift.x(), -shift.y());
        // Remove points outside snug bbox (simple pre-clip for performance).
        pts.erase(
            std::remove_if(pts.begin(), pts.end(),
                           [&](const Point &p) { return !snug_shifted.contains(p); }),
            pts.end());
    }

    if (pts.size() < 2)
        return;

    Polyline polyline;
    polyline.points = std::move(pts);

    Polylines polylines = intersection_pl(polyline, expolygon);

    if (!polylines.empty()) {
        Polylines chained;
        if (params.dont_connect() || params.density > 0.5 || polylines.size() <= 1)
            chained = chain_polylines(std::move(polylines));
        else
            connect_infill(std::move(polylines), expolygon, chained, this->spacing, params);

        for (Polyline &pl : chained) {
            pl.translate(shift.x(), shift.y());
            pl.rotate(direction.first);
        }
        append(polylines_out, std::move(chained));
    }
}

} // namespace Slic3r
