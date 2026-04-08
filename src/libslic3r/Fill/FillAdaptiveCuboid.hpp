// Adaptive Cuboid infill — a variant of Adaptive Cubic with anisotropic (cuboid) cells.
// Cells are stretched along the Z axis by a configurable ratio, producing taller cells
// that use less material vertically while maintaining horizontal structural strength.
// Reuses the FillAdaptive octree infrastructure with a Z-scaling transform.

#ifndef slic3r_FillAdaptiveCuboid_hpp_
#define slic3r_FillAdaptiveCuboid_hpp_

#include "FillAdaptive.hpp"

namespace Slic3r {
namespace FillAdaptiveCuboid {

class Filler : public FillAdaptive::Filler
{
public:
    ~Filler() override {}

protected:
    Fill* clone() const override { return new Filler(*this); }
    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;
};

} // namespace FillAdaptiveCuboid
} // namespace Slic3r

#endif // slic3r_FillAdaptiveCuboid_hpp_
