// Adaptive Cuboid infill implementation.
// Delegates to FillAdaptive::Filler with a Z-scaling transform applied:
// - Before generating infill lines, the Z position is divided by z_stretch,
//   making the octree produce lines as if the object were shorter in Z.
// - This effectively stretches octree cells along Z, creating cuboid (non-cube) cells.

#include "FillAdaptiveCuboid.hpp"

namespace Slic3r {
namespace FillAdaptiveCuboid {

void Filler::_fill_surface_single(
    const FillParams              &params,
    unsigned int                   thickness_layers,
    const std::pair<float, Point> &direction,
    ExPolygon                      expolygon,
    Polylines                     &polylines_out)
{
    // Apply anisotropic Z scaling: divide the Z position by the stretch factor
    // so the octree sees a compressed Z coordinate, effectively stretching cells vertically.
    const coordf_t original_z = this->z;
    if (this->z_stretch > 1.0)
        this->z = original_z / this->z_stretch;

    // Delegate to the parent FillAdaptive::Filler implementation
    FillAdaptive::Filler::_fill_surface_single(params, thickness_layers, direction, std::move(expolygon), polylines_out);

    // Restore original Z
    this->z = original_z;
}

} // namespace FillAdaptiveCuboid
} // namespace Slic3r
