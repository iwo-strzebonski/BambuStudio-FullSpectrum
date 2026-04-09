#ifndef slic3r_FillArchimedeanSpiral3D_hpp_
#define slic3r_FillArchimedeanSpiral3D_hpp_

#include "../libslic3r.h"
#include "FillPlanePath.hpp"

namespace Slic3r {

// Archimedean Spiral 3D infill.
//
// Generates a single outward Archimedean spiral per layer, alternating the
// winding direction (CW / CCW) every layer.  The spiral start angle is offset
// by the golden angle (≈ 137.5°) multiplied by the layer index, so consecutive
// layers begin at a different angular position.  The resulting cross-layer
// path intersections form a 3-D interlocking weave that improves inter-layer
// adhesion without requiring retractions within any single layer.
//
// Properties:
//   - Near-zero in-layer retractions (single continuous path per layer)
//   - Good inter-layer bonding (alternating direction cross-hatch)
//   - Poor bridging (no straight unsupported segments)
//   - Good for flexible materials
class FillArchimedeanSpiral3D : public FillPlanePath
{
public:
    Fill* clone() const override { return new FillArchimedeanSpiral3D(*this); }
    ~FillArchimedeanSpiral3D() override = default;

protected:
    bool centered() const override { return true; }

    void generate(coord_t min_x, coord_t min_y, coord_t max_x, coord_t max_y,
                  const double resolution, InfillPolylineOutput &output) override;
};

} // namespace Slic3r

#endif // slic3r_FillArchimedeanSpiral3D_hpp_
