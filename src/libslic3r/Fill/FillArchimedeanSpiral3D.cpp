// FillArchimedeanSpiral3D.cpp
//
// Archimedean Spiral 3D infill for BambuStudio-FullSpectrum.
//
// Each layer contains a single continuous Archimedean spiral (r = a + b*θ)
// that starts at the centre and unwinds to the bounding-circle radius.
//
// Two 3-D effects are layered on top of the basic 2-D spiral:
//
//   1. Direction alternation — the winding sense (CCW/CW) flips every layer.
//      The spiral on layer N crosses the spiral on layer N+1 at every ring,
//      creating Z-direction interlocks analogous to those in FillGyroid but
//      with a much simpler, single-path topology.
//
//   2. Golden-angle phase offset — the angular start of each spiral is
//      advanced by the golden angle (π·(3−√5) ≈ 137.508°) times the layer
//      index.  The golden angle maximises angular separation between
//      successive layers, so the outermost (densest) ring of one layer is
//      never angularly co-located with the same ring of the adjacent layer.
//      This smears the small gap at r=0 around the full circumference over
//      several layers, preventing a persistent thin column of under-filled
//      material at the centre.

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include "FillArchimedeanSpiral3D.hpp"

#include <cmath>
#include <algorithm>

namespace Slic3r {

// Generate an Archimedean spiral (r = a + b*θ) into an InfillPolylineOutput.
//
// InfillPolylineClipper is defined only in FillPlanePath.cpp (it is not
// exposed in the header), so we cannot use the template + static_cast trick
// employed by FillArchimedeanChords.  Instead we call add_point() through the
// base-class interface.  The intermediate bbox-clipping optimisation is
// therefore skipped, but correctness is fully preserved: FillPlanePath::
// _fill_surface_single() calls intersection_pl() after generate(), which clips
// the result to the actual expolygon boundary regardless.
static void generate_archimedean_spiral_3d(
    coord_t               min_x,
    coord_t               min_y,
    coord_t               max_x,
    coord_t               max_y,
    const double          resolution,
    const bool            reverse,
    const double          phase,
    InfillPolylineOutput &output)
{
    // Radius required to fully cover the bounding rectangle (diagonal + margin).
    const coordf_t rmax = std::sqrt(
        coordf_t(max_x) * coordf_t(max_x) +
        coordf_t(max_y) * coordf_t(max_y)) * std::sqrt(2.) + 1.5;

    // Spiral parameters: r = a + b*theta  (a=1 offsets centre gap, b=1/2π
    // makes adjacent rings exactly 1 unit apart in normalised coordinates).
    const coordf_t a   = 1.;
    const coordf_t b   = 1. / (2. * M_PI);
    const double   dir = reverse ? -1.0 : 1.0;

    coordf_t theta = 0.;
    coordf_t r     = a; // matches the second seed point below

    // Seed: origin anchor, then the first spiral point at r=a, phase-rotated.
    output.add_point({ 0., 0. });
    output.add_point({ std::cos(phase), std::sin(phase) });

    while (r < rmax) {
        // Advance θ by the angle that produces a chord length equal to
        // 'resolution' at the current radius (same adaptive discretisation
        // as FillArchimedeanChords).
        const double acos_arg = std::max(-1.0, std::min(1.0, 1.0 - resolution / r));
        theta += 2.0 * std::acos(acos_arg);
        r      = a + b * theta;

        const double eff_theta = phase + dir * theta;
        output.add_point({ r * std::cos(eff_theta), r * std::sin(eff_theta) });
    }
}

void FillArchimedeanSpiral3D::generate(
    coord_t min_x, coord_t min_y, coord_t max_x, coord_t max_y,
    const double resolution, InfillPolylineOutput &output)
{
    // Odd layers get the reversed (CW) winding direction, creating an
    // alternating CCW/CW interlocking structure across the Z axis.
    const bool   reverse = (this->layer_id % 2) != 0;

    // Golden-angle phase offset: π·(3−√5) ≈ 2.39996 rad ≈ 137.508°.
    // Successive layers start at maximally-spread angular positions.
    const double phase = static_cast<double>(this->layer_id) * M_PI * (3.0 - std::sqrt(5.0));

    generate_archimedean_spiral_3d(min_x, min_y, max_x, max_y, resolution,
                                   reverse, phase, output);
}

} // namespace Slic3r
