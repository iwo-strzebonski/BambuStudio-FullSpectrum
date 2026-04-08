#ifndef slic3r_FillSchwartzP_hpp_
#define slic3r_FillSchwartzP_hpp_

#include "../libslic3r.h"

#include "FillBase.hpp"

namespace Slic3r {

// Schwartz P TPMS infill.  The minimal surface is defined by cos(x)+cos(y)+cos(z)=0,
// which produces a periodic bicontinuous foam ideal for structural infill.
class FillSchwartzP : public Fill
{
public:
    FillSchwartzP() {}
    Fill* clone() const override { return new FillSchwartzP(*this); }

    bool use_bridge_flow() const override { return false; }
    bool is_self_crossing() override { return false; }

    // Correction applied to infill angle (degrees)
    static constexpr float CorrectionAngle = -45.;

    // Density adjustment for a good weight percentage.
    static constexpr double DensityAdjust = 2.0;

    // Upper resolution tolerance (mm^-2)
    static constexpr double PatternTolerance = 0.2;

protected:
    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillSchwartzP_hpp_
