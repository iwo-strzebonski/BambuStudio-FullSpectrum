#ifndef slic3r_FillHilbertCurve3D_hpp_
#define slic3r_FillHilbertCurve3D_hpp_

#include "../libslic3r.h"
#include "FillBase.hpp"

namespace Slic3r {

// 3D Hilbert Curve infill.
//
// Computes the XY cross-section of a true 3D Hilbert space-filling curve for
// each layer.  Because the 3D curve's recursion state machine reorients the
// sub-curve differently at every Z-level, adjacent layers receive complementary
// path orientations — producing Z-crossing interlocks that dramatically improve
// inter-layer bonding and impact resistance compared with the 2D Hilbert curve
// (which is simply re-used unchanged on every layer).
//
// Implementation is completely stateless: the Z coordinate is derived from
// Fill::layer_id so no cross-layer shared state is required.
//
// The pattern repeats every 2^order layers; the order is chosen to be the
// smallest power-of-two cube side that covers both the XY bounding box and a
// configurable maximum-order cap (MAX_ORDER = 6, i.e. 64×64×64 = 262 144 curve
// nodes per repetition period).
class FillHilbertCurve3D : public Fill
{
public:
    FillHilbertCurve3D() = default;
    Fill* clone() const override { return new FillHilbertCurve3D(*this); }
    ~FillHilbertCurve3D() override = default;

    bool is_self_crossing() override { return false; }
    bool has_consistent_pattern() const override { return true; }

    // Maximum curve order.  Order N covers a 2^N × 2^N × 2^N cube → caps the
    // O(sz³) iteration at 262 144 nodes per layer call.
    static constexpr int MAX_ORDER = 6;

protected:
    float _layer_angle(size_t /*idx*/) const override { return 0.f; }

    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillHilbertCurve3D_hpp_
