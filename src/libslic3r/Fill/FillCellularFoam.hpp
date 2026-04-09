#ifndef slic3r_FillCellularFoam_hpp_
#define slic3r_FillCellularFoam_hpp_

#include "../libslic3r.h"
#include "FillBase.hpp"

namespace Slic3r {

// Cellular Foam infill.
//
// Deliberately randomized cell sizes and positions within a density budget.
// Mimics closed-cell foam microstructure: a Voronoi tessellation seeded by a
// stochastic multi-scale point distribution produces irregular polygon cells
// of varying sizes — exactly like a cross-section through rigid closed-cell
// foam (PU foam, EVA foam, etc.).
//
// Why it damps vibrations better than periodic patterns:
//   Periodic patterns (grid, honeycomb, gyroid) have a dominant spatial
//   frequency equal to their cell pitch.  Any mechanical excitation close to
//   the resonant frequency of that pitch propagates through the part nearly
//   without attenuation.  The intentional size/position disorder here breaks
//   that periodicity; energy is scattered across a continuum of spatial
//   frequencies rather than channelled into a single resonant mode.
//
// Algorithm overview:
//   1. Seed a PCG32 PRNG from (layer_id ^ bbox_hash) for fully reproducible
//      but layer-unique random patterns.
//   2. Generate multi-scale seeds in two passes over the bounding box:
//      - Coarse pass  (step ≈ 2.0 × base_spacing): accepted with probability
//        CoarseAcceptProb — produces the larger "bubble" cells.
//      - Fine pass    (step ≈ 0.85 × base_spacing): accepted with probability
//        FineAcceptProb — fills gaps between coarse bubbles with smaller cells.
//      Both passes apply random jitter to the grid position so no two layers
//      share the same cell layout.  A spatial hash enforces a minimum
//      seed-to-seed separation so cells never collapse to zero width.
//   3. Build a 2-D Voronoi diagram from the accepted seeds (boost::polygon).
//   4. Clip all Voronoi edges (including infinite ones) to the ExPolygon.
//   5. Discard degenerate segments shorter than MinSegFraction × base_spacing.
//   6. Sort the remaining polylines with chain_polylines() and emit them.
//
// Density control:
//   base_spacing = this->spacing / params.density
//   Smaller density → larger base_spacing → fewer seeds → fewer cell walls →
//   less material, consistent with how every other infill pattern works.
class FillCellularFoam : public Fill
{
public:
    FillCellularFoam() {}
    Fill* clone() const override { return new FillCellularFoam(*this); }

    bool use_bridge_flow() const override { return false; }
    bool is_self_crossing()        override { return false; }

    // Grid step for the large-cell (coarse) pass as a multiple of base_spacing.
    static constexpr double CoarseScale       = 2.0;

    // Grid step for the small-cell (fine) pass as a multiple of base_spacing.
    static constexpr double FineScale         = 0.85;

    // Fraction of coarse grid candidates that produce a seed.
    // ~55 % acceptance at the coarse level gives well-separated large cells.
    static constexpr double CoarseAcceptProb  = 0.55;

    // Fraction of fine grid candidates that produce a seed.
    // ~38 % acceptance fills the gaps without over-crowding.
    static constexpr double FineAcceptProb    = 0.38;

    // Maximum jitter as a fraction of the respective grid step.
    // Coarse cells are displaced by up to ± JitterFrac × coarse_step in X and Y.
    static constexpr double JitterFrac        = 0.42;

    // Minimum seed-to-seed distance as a fraction of base_spacing.
    // Prevents degenerate zero-width cell walls.
    static constexpr double MinSepFrac        = 0.45;

    // Voronoi edges shorter than this fraction of base_spacing are discarded.
    static constexpr double MinSegFrac        = 0.25;

protected:
    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillCellularFoam_hpp_
