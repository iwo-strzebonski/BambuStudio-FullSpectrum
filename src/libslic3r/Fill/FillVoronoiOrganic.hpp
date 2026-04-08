#ifndef slic3r_FillVoronoiOrganic_hpp_
#define slic3r_FillVoronoiOrganic_hpp_

#include "../libslic3r.h"
#include "FillBase.hpp"

namespace Slic3r {

// Voronoi Organic infill.
//
// Voronoi tessellation seeded by a stress-following point distribution.
// Seeds are placed at high density near the shell boundaries (where bending
// stress concentrations occur) and at lower density in the bulk interior —
// mimicking the principal-stress-driven material distribution used in
// aerospace composite lattices.
//
// Algorithm overview:
//   1. Generate seed points using a variable-spacing Poisson-disk process
//      driven by the squared distance to the nearest polygon edge.
//   2. Build a 2-D Voronoi diagram from those seeds (boost::polygon).
//   3. Clip all Voronoi edges (including infinite ones) to the ExPolygon.
//   4. Sort the resulting segments with chain_polylines() for print
//      efficiency and emit them as Polylines.
//
// The z coordinate is mixed into the seed jitter so that the organic-looking
// cell pattern shifts gradually between layers.
class FillVoronoiOrganic : public Fill
{
public:
    FillVoronoiOrganic() {}
    Fill* clone() const override { return new FillVoronoiOrganic(*this); }

    bool use_bridge_flow() const override { return false; }
    bool is_self_crossing()        override { return false; }

    // Seed spacing near the wall as a fraction of base_spacing.
    // Smaller  → denser cells at the stress concentration zone.
    static constexpr double DenseFactor  = 0.45;

    // Seed spacing in the stress-free interior as a fraction of base_spacing.
    static constexpr double SparseFactor = 2.2;

    // Width of the transition zone expressed as a multiple of base_spacing.
    // At distance > FalloffFactor*base_spacing the density is fully sparse.
    static constexpr double FalloffFactor = 2.0;

    // Maximum jitter amplitude relative to dense_spacing.
    // Keeps the pattern organic without disrupting the density gradient.
    static constexpr double JitterAmplitude = 0.4;

protected:
    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillVoronoiOrganic_hpp_
