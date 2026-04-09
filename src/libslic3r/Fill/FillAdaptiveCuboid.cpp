// Adaptive Cuboid infill implementation.
// The cuboid Z-stretch is applied in FillAdaptive::Filler::_fill_surface_single
// via the z_stretch field on FillContext. This subclass exists only as a
// distinct type for enum dispatch; z_stretch is set from adaptive_cuboid_z_ratio
// config by Fill.cpp before calling fill_surface.

#include "FillAdaptiveCuboid.hpp"

namespace Slic3r {
namespace FillAdaptiveCuboid {

// No overrides needed. FillAdaptive::Filler handles z_stretch generically.

} // namespace FillAdaptiveCuboid
} // namespace Slic3r
