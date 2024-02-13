#include <util/vector.h>

namespace lightfold {

	class Point {
    public:
        static const int dim = 3;

        // Point Public Methods
        Point() : x(0), y(0), z(0) {}
        Point(Float x, Float y, Float z) : x(x), y(y), z(z) {}

        bool operator==(Point c) const {
            return x == c.x && y == c.y && z == c.z;
        }
        bool operator!=(Point c) const {
            return x != c.x || y != c.y || z != c.z;
        }

        Float operator[](int i) const {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
        Float& operator[](int i) {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        // Point Public Members
        Float x, y, z;
	};

    Point transport(Point p, Vector v, Float multiplier = 1) {
        return Point(p.x + multiplier * v.x, p.y + multiplier * v.y, p.z + multiplier * v.z);
    }

    Vector fromTo(Point p, Point q) {
        return Vector(q.x - p.x, q.y - p.y, q.z - p.z);
    }

} // namespace lightfold