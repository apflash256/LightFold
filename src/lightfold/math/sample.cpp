#include <math/sample.h>

namespace lightfold {

    Point2f ConcentricSampleDisk(const Point2f& u) {
        Point2f uOffset = 2.f * u - GeoVector2f(1, 1);
        if (uOffset.x == 0 && uOffset.y == 0)
            return Point2f(0, 0);
        float theta, r;
        if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
            r = uOffset.x;
            theta = PiOver4 * (uOffset.y / uOffset.x);
        }
        else {
            r = uOffset.y;
            theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
        }
        return r * Point2f(std::cos(theta), std::sin(theta));
    }

} // namespace lightfold