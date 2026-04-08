// Schwartz P TPMS infill.
// The surface is defined by cos(x) + cos(y) + cos(z) = 0.
// The wave function reduces to: y = arccos(clamp(-cos(z) - cos(x), -1, 1))
// which is derived by solving for y in the TPMS equation at a given (x, z).
// Implementation follows the same wave-generation pipeline as FillGyroid.

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include <cmath>
#include <algorithm>

#include "FillSchwartzP.hpp"

namespace Slic3r {

// Schwartz P surface: cos(x) + cos(y) + cos(z) = 0
// Solved for y: y = arccos(clamp(-cos(z) - cos(x), -1, 1))
// z_cos = cos(z) is pre-computed for performance.
// The `flip` argument mirrors the solution to the negative branch (y -> 2*PI - y).
static inline double f_schwartz(double x, double z_cos, bool flip)
{
    double val = std::clamp(-z_cos - std::cos(x), -1.0, 1.0);
    double y   = std::acos(val);
    return flip ? (2.0 * M_PI - y) : y;
}

static inline Polyline make_wave_schwartz(
    const std::vector<Vec2d>& one_period, double width, double height, double offset,
    double scaleFactor, double z_cos, bool flip)
{
    std::vector<Vec2d> points = one_period;
    double period = points.back()(0);
    if (width != period) {
        points.reserve(one_period.size() * size_t(std::floor(width / period)));
        points.pop_back();

        size_t n = points.size();
        do {
            points.emplace_back(points[points.size() - n].x() + period,
                                points[points.size() - n].y());
        } while (points.back()(0) < width - EPSILON);

        points.emplace_back(Vec2d(width, f_schwartz(width, z_cos, flip)));
    }

    Polyline polyline;
    polyline.points.reserve(points.size());
    for (auto& point : points) {
        point(1) += offset;
        point(1) = std::clamp(double(point.y()), 0., height);
        polyline.points.emplace_back((point * scaleFactor).cast<coord_t>());
    }
    return polyline;
}

static std::vector<Vec2d> make_one_period_schwartz(
    double width, double scaleFactor, double z_cos, bool flip, double tolerance)
{
    std::vector<Vec2d> points;
    double dx    = M_PI_2;
    double limit = std::min(2.0 * M_PI, width);
    points.reserve(coord_t(std::ceil(limit / tolerance / 3)));

    for (double x = 0.; x < limit - EPSILON; x += dx)
        points.emplace_back(Vec2d(x, f_schwartz(x, z_cos, flip)));
    points.emplace_back(Vec2d(limit, f_schwartz(limit, z_cos, flip)));

    // Adaptive refinement up to requested tolerance.
    for (;;) {
        size_t size = points.size();
        for (unsigned int i = 1; i < size; ++i) {
            auto& lp = points[i - 1];
            auto& rp = points[i];
            double x  = lp(0) + (rp(0) - lp(0)) / 2.0;
            double y  = f_schwartz(x, z_cos, flip);
            Vec2d  ip = {x, y};
            if (std::abs(cross2(Vec2d(ip - lp), Vec2d(ip - rp))) > sqr(tolerance))
                points.emplace_back(std::move(ip));
        }
        if (size == points.size())
            break;
        std::sort(points.begin(), points.end(),
                  [](const Vec2d &a, const Vec2d &b) { return a(0) < b(0); });
    }

    return points;
}

static Polylines make_schwartz_waves(double gridZ, double density_adjusted, double line_spacing, double width, double height)
{
    const double scaleFactor = scale_(line_spacing) / density_adjusted;
    const double tolerance   = std::min(line_spacing / 2.0, FillSchwartzP::PatternTolerance) / unscale<double>(scaleFactor);

    const double z     = gridZ / scaleFactor;
    const double z_cos = std::cos(z);

    // Odd rows: y = arccos(-cos(z)-cos(x)), even rows: mirrored branch
    std::vector<Vec2d> one_period_odd  = make_one_period_schwartz(width, scaleFactor, z_cos, false, tolerance);
    std::vector<Vec2d> one_period_even = make_one_period_schwartz(width, scaleFactor, z_cos, true,  tolerance);

    Polylines result;
    for (double y0 = 0.; y0 < height + EPSILON; y0 += M_PI) {
        result.emplace_back(make_wave_schwartz(one_period_odd,  width, height, y0, scaleFactor, z_cos, false));
        y0 += M_PI;
        if (y0 < height + EPSILON)
            result.emplace_back(make_wave_schwartz(one_period_even, width, height, y0, scaleFactor, z_cos, true));
    }

    return result;
}

// FIXME: needed to fix build on Mac buildserver
constexpr double FillSchwartzP::PatternTolerance;

void FillSchwartzP::_fill_surface_single(
    const FillParams                &params,
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction,
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    auto infill_angle = float(this->angle + (CorrectionAngle * 2 * M_PI) / 360.);
    if (std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    BoundingBox bb = expolygon.contour.bounding_box();
    double      density_adjusted = std::max(0., params.density * DensityAdjust / params.multiline);
    coord_t     distance         = coord_t(scale_(this->spacing) / density_adjusted);

    bb.merge(align_to_grid(bb.min, Point(2 * M_PI * distance, 2 * M_PI * distance)));

    Polylines polylines = make_schwartz_waves(
        scale_(this->z),
        density_adjusted,
        this->spacing,
        std::ceil(bb.size()(0) / distance) + 1.,
        std::ceil(bb.size()(1) / distance) + 1.);

    for (Polyline &pl : polylines)
        pl.translate(bb.min);

    multiline_fill(polylines, params, spacing);

    polylines = intersection_pl(polylines, expolygon);

    if (! polylines.empty()) {
        const double minlength = scale_(0.8 * this->spacing);
        polylines.erase(
            std::remove_if(polylines.begin(), polylines.end(),
                           [minlength](const Polyline &pl) { return pl.length() < minlength; }),
            polylines.end());
    }

    if (! polylines.empty()) {
        size_t polylines_out_first_idx = polylines_out.size();
        if (params.dont_connect())
            append(polylines_out, chain_polylines(polylines));
        else
            this->connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);

        if (std::abs(infill_angle) >= EPSILON) {
            for (auto it = polylines_out.begin() + polylines_out_first_idx; it != polylines_out.end(); ++it)
                it->rotate(infill_angle);
        }
    }
}

} // namespace Slic3r
