#pragma once

#include <math/vecmath.h>

namespace lightfold {

    class Filter {
    public:
        // Filter Public Methods
        virtual ~Filter() {}
        Filter(const GeoVector2f& radius)
            : radius(radius), invRadius(GeoVector2f(1 / radius.x, 1 / radius.y)) { }
        virtual float Evaluate(const Point2f& p) const = 0;

        // Filter Public Data
        const GeoVector2f radius, invRadius;
    };

    class MitchellFilter : public Filter {
    public:
        // MitchellFilter Public Methods
        MitchellFilter(const GeoVector2f& radius, float B, float C)
            : Filter(radius), B(B), C(C) { }
        float Evaluate(const Point2f& p) const;
        float Mitchell1D(float x) const {
            x = std::abs(2 * x);
            if (x > 1)
                return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x + (-12 * B - 48 * C) * x + (8 * B + 24 * C)) * (1.f / 6.f);
            else
                return ((12 - 9 * B - 6 * C) * x * x * x + (-18 + 12 * B + 6 * C) * x * x + (6 - 2 * B)) * (1.f / 6.f);
        }

    private:
        // MitchellFilter Private Data
        const float B, C;
    };

} // namespace lightfold