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
    // _fill_surface_single is inherited from FillAdaptive::Filler.
    // The cuboid Z-stretch is handled there via this->z_stretch,
    // which is set from adaptive_cuboid_z_ratio config in Fill.cpp.
};

} // namespace FillAdaptiveCuboid
} // namespace Slic3r

#endif // slic3r_FillAdaptiveCuboid_hpp_
