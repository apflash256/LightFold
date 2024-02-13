#pragma once

#include <util/math.h>

namespace lightfold {

    template <typename T>
    class Point2 {
    public:
        static const int dim = 2;

        // Point2 Public Methods
        Point2() : x(T()), y(T()) {}
        Point2(T x, T y) : x(x), y(y) {}

        template <typename U>
        operator Point2<U>() {
            return Point2<U>(U(x), U(y));
        }

        bool operator==(Point2<T> c) const {
            return x == c.x && y == c.y;
        }
        bool operator!=(Point2<T> c) const {
            return x != c.x || y != c.y;
        }

        T operator[](int i) const {
            return (i == 0) ? x : y;
        }
        T& operator[](int i) {
            return (i == 0) ? x : y;
        }

        // Point2 Public Members
        T x, y;
    };

    template <typename T>
    class Vector2 {
    public:
        static const int dim = 2;

        // Vector2 Public Methods
        Vector2() : x(T()), y(T()) {}
        Vector2(T x, T y) : x(x), y(y) {}

        template <typename U>
        operator Vector2<U>() {
            return Vector2<U>(U(x), U(y));
        }

        bool operator==(Vector2<T> c) const {
            return x == c.x && y == c.y;
        }
        bool operator!=(Vector2<T> c) const {
            return x != c.x || y != c.y;
        }

        T operator[](int i) const {
            return (i == 0) ? x : y;
        }
        T& operator[](int i) {
            return (i == 0) ? x : y;
        }

        // Vector2 Public Members
        T x, y;
    };

    class Point;
    class Vector;

    Point transport(Point p, Vector v, Float multiplier = 1);
    Vector fromTo(Point p, Point q);

    class Vector {
    public:
        static const int dim = 3;

        // Vector Public Methods
        Vector() : x(0), y(0), z(0) {}
        Vector(Float x, Float y, Float z) : x(x), y(y), z(z) {}

        bool operator==(Vector c) const {
            return x == c.x && y == c.y && z == c.z;
        }
        bool operator!=(Vector c) const {
            return x != c.x || y != c.y || z != c.z;
        }

        Float operator[](int i) const {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
        Float &operator[](int i) {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        // Vector Public Members
        Float x, y, z;
    };

    class Normal {
    public:
        static const int dim = 3;

        // Normal Public Methods
        Normal() : x(0), y(0), z(0) {}
        Normal(Float x, Float y, Float z) : x(x), y(y), z(z) {}

        bool operator==(Normal c) const {
            return x == c.x && y == c.y && z == c.z;
        }
        bool operator!=(Normal c) const {
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

        // Normal Public Members
        Float x, y, z;
    };

    class Ray {
    public:
        // Ray Public Methods
        Ray() : origin(Point()), dir(Vector()) {}
        Ray(Point origin, Vector dir) : origin(origin), dir(dir) {}

        // Ray Public Members
        Point origin;
        Vector dir;
    };

    class Transformation {
    public:
        Transformation() : matrix{ {1,0,0},{0,1,0},{0,0,1} }, invmat{ {1,0,0},{0,1,0},{0,0,1} } {}

        Transformation(const Float m[3][3]) {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    matrix[i][j] = m[i][j];
            invmat[0][0] = m[1][1] * m[2][2] - m[2][1] * m[1][2];
            invmat[0][1] = m[2][1] * m[0][2] - m[0][1] * m[2][2];
            invmat[0][2] = m[0][1] * m[1][2] - m[1][1] * m[0][2];
            Float invdet = 1.0 /
                (invmat[0][0] * m[0][0] + invmat[0][1] * m[1][0] + invmat[0][2] * m[2][0]);
            invmat[0][0] *= invdet;
            invmat[0][1] *= invdet;
            invmat[0][2] *= invdet;
            invmat[1][0] = (m[1][2] * m[2][0] - m[2][2] * m[1][0]) * invdet;
            invmat[1][1] = (m[2][2] * m[0][0] - m[0][2] * m[2][0]) * invdet;
            invmat[1][2] = (m[0][2] * m[1][0] - m[1][2] * m[0][0]) * invdet;
            invmat[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
            invmat[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
            invmat[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
        }

        Transformation(const Float m[3][3], const Float minv[3][3]) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    matrix[i][j] = m[i][j];
                    invmat[i][j] = minv[i][j];
                }
            }
        }

        // Transformation Public Members
        Float matrix[3][3];
        Float invmat[3][3];
    };

} // namespace lightfold