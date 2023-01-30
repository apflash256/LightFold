#pragma once

#include <core/floating.h>
#include <math/math.h>

#include <array>
#include <type_traits>
#include <algorithm>

namespace lightfold {

    template <typename T>
    class Vector2;
    template <typename T>
    class Vector3;
    using Vector2f = Vector2<float>;
    using Vector2i = Vector2<int>;
    using Vector3f = Vector3<float>;
    using Vector3i = Vector3<int>;

    template <typename T>
    class Point2;
    template <typename T>
    class Point3;
    using Point2f = Point2<float>;
    using Point2i = Point2<int>;
    using Point3f = Point3<float>;
    using Point3i = Point3<int>;
    class Point3fi;

    template <typename T>
    class Tangent3;
    using Tangent3f = Tangent3<float>;
    class Tangent3fi;

    template <typename T>
    class Normal3;
    using Normal3f = Normal3<float>;

    template <typename T>
    class Bounds2;
    template <typename T>
    class Bounds3;
    using Bounds2f = Bounds2<float>;
    using Bounds2i = Bounds2<int>;
    using Bounds3f = Bounds3<float>;
    using Bounds3i = Bounds3<int>;

    namespace {
        // TupleLength Definition
        template <typename T>
        struct TupleLength {
            using type = float; // use float for int, long, ...
        };

        template <>
        struct TupleLength<double> {
            using type = double;
        };

        template <>
        struct TupleLength<long double> {
            using type = long double;
        };
    } // anonymous namespace

    // Class Definitions
    template <typename T>
    class Tuple2 {
    public:
        static const int nDimensions = 2;

        // Tuple2 Public Methods
        Tuple2() = default;
        Tuple2(T x, T y) : x(x), y(y) {}

        bool operator==(Tuple2<T> c) const {
            return x == c.x && y == c.y;
        }
        bool operator!=(Tuple2<T> c) const {
            return x != c.x || y != c.y;
        }

        T operator[](int i) const {
            return (i == 0) ? x : y;
        }
        T& operator[](int i) {
            return (i == 0) ? x : y;
        }

        // Tuple2 Public Members
        T x, y;
    };

    template <typename T>
    class Tuple3 {
    public:
        static const int nDimensions = 3;

        // Tuple3 Public Methods
        Tuple3() = default;
        Tuple3(T x, T y, T z) : x(x), y(y), z(z) { }

        bool operator==(Tuple3<T> c) const {
            return x == c.x && y == c.y && z == c.z;
        }
        bool operator!=(Tuple3<T> c) const {
            return x != c.x || y != c.y || z != c.z;
        }

        T operator[](int i) const {
            if (i == 0)
                return x;
            if (i == 1)
                return y;
            return z;
        }
        T& operator[](int i) {
            if (i == 0)
                return x;
            if (i == 1)
                return y;
            return z;
        }

        // Tuple3 Public Members
        T x, y, z;
    };

    template <typename T>
    class Point2 : public Tuple2<T> {
    public:
        // Point2 Public Methods
        Point2() = default;
        Point2(T x, T y) : Tuple2<T>(x, y) {}

        template <typename U>
        explicit Point2(Tuple2<U> p) : Tuple2<T>(T(p.x), T(p.y)) {}

        // We cannot add two points, but we can add a point and a vector.
        // Here we specify the vectors to be tangent vectors.
        template <typename U>
        auto operator+(Vector2<U> v) const->Point2<decltype(T{} + U{}) > {
            return { x + v.x, y + v.y };
        }
        template <typename U>
        Point2<T>& operator+=(Vector2<U> v) {
            x += v.x;
            y += v.y;
            return *this;
        }
        template <typename U>
        auto operator-(Point2<U> p) const->Vector2<decltype(T{} - U{}) > {
            return { x - p.x, y - p.y };
        }
        template <typename U>
        auto operator-(Vector2<U> v) const->Point2<decltype(T{} - U{}) > {
            return { x - v.x, y - v.y };
        }
        template <typename U>
        Point2<T>& operator-=(Vector2<U> v) {
            x -= v.x;
            y -= v.y;
            return *this;
        }

        using Tuple2<T>::x;
        using Tuple2<T>::y;
    };

    template <typename T>
    class Point3 : public Tuple3<T> {
    public:
        // Point3 Public Methods
        Point3() = default;
        Point3(T x, T y, T z) : Tuple3<T>(x, y, z) {}

        template <typename U>
        explicit Point3(Tuple3<U> p) : Tuple3<T>(T(p.x), T(p.y), T(p.z)) {}

        // We cannot add two points, but we can add a point and a vector.
        // Here we specify the vectors to be tangent vectors.
        template <typename U>
        auto operator+(Tangent3<U> v) const->Point3<decltype(T{} + U{}) > {
            return { x + v.x, y + v.y, z + v.z };
        }
        template <typename U>
        Point3<T>& operator+=(Tangent3<U> v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        template <typename U>
        auto operator-(Point3<U> p) const->Tangent3<decltype(T{} - U{}) > {
            return { x - p.x, y - p.y, z - p.z };
        }
        template <typename U>
        auto operator-(Tangent3<U> v) const->Point3<decltype(T{} - U{}) > {
            return { x - v.x, y - v.y, z - v.z };
        }
        template <typename U>
        Point3<T>& operator-=(Tangent3<U> v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        using Tuple3<T>::x;
        using Tuple3<T>::y;
        using Tuple3<T>::z;
    };

    template <typename T>
    class Vector2 : public Tuple2<T> {
    public:
        // Vector2 Public Methods
        Vector2() = default;
        Vector2(T x, T y) : Tuple2<T>(x, y) {}

        template <typename U>
        auto operator+(Vector2<U> c) const->Vector2<decltype(T{} + U{}) > {
            return { x + c.x, y + c.y };
        }
        template <typename U>
        Vector2<T>& operator+=(Vector2<U> c) {
            x += c.x;
            y += c.y;
            return static_cast<Vector2<T> &>(*this);
        }
        template <typename U>
        auto operator-(Vector2<U> c) const->Vector2<decltype(T{} - U{}) > {
            return { x - c.x, y - c.y };
        }
        template <typename U>
        Vector2<T>& operator-=(Vector2<U> c) {
            x -= c.x;
            y -= c.y;
            return static_cast<Vector2<T> &>(*this);
        }
        template <typename U>
        Vector2<T>& operator*=(U s) {
            x *= s;
            y *= s;
            return static_cast<Vector2<T> &>(*this);
        }
        template <typename U>
        Vector2<T>& operator/=(U d) {
            x /= d;
            y /= d;
            return static_cast<Vector2<T> &>(*this);
        }

        template <typename U>
        explicit Vector2(Vector2<U> v) : Tuple2<T>(T(v.x), T(v.y)) {}
        template <typename U>
        explicit Vector2(Point2<U> p) : Tuple2<T>(T(p.x), T(p.y)) {}

        using Tuple2<T>::x;
        using Tuple2<T>::y;
    };

    template <typename T>
    class Vector3 : public Tuple3<T> {
    public:
        // Vector3 Public Methods
        Vector3() = default;
        Vector3(T x, T y, T z) : Tuple3<T>(x, y, z) {}

        template <typename U>
        auto operator+(Vector3<U> c) const->Vector3<decltype(T{} + U{}) > {
            return { x + c.x, y + c.y, z + c.z };
        }
        template <typename U>
        Vector3<T>& operator+=(Vector3<U> c) {
            x += c.x;
            y += c.y;
            z += c.z;
            return static_cast<Vector3<T> &>(*this);
        }
        template <typename U>
        auto operator-(Vector3<U> c) const->Vector3<decltype(T{} - U{}) > {
            return { x - c.x, y - c.y, z - c.z };
        }
        template <typename U>
        Vector3<T>& operator-=(Vector3<U> c) {
            x -= c.x;
            y -= c.y;
            z -= c.z;
            return static_cast<Vector3<T> &>(*this);
        }
        template <typename U>
        Vector3<T>& operator*=(U s) {
            x *= s;
            y *= s;
            z *= s;
            return static_cast<Vector3<T> &>(*this);
        }
        template <typename U>
        Vector3<T>& operator/=(U d) {
            x /= d;
            y /= d;
            z /= d;
            return static_cast<Vector3<T> &>(*this);
        }

        template <typename U>
        explicit Vector3(Vector3<U> v) : Tuple3<T>(T(v.x), T(v.y), T(v.z)) {}
        template <typename U>
        explicit Vector3(Point3<U> p) : Tuple3<T>(T(p.x), T(p.y), T(p.z)) {}
        template <typename U>
        explicit Vector3(Tangent3<U> n) : Tuple3<T>(T(n.x), T(n.y), T(n.z)) {}
        template <typename U>
        explicit Vector3(Normal3<U> n) : Tuple3<T>(T(n.x), T(n.y), T(n.z)) {}

        using Tuple3<T>::x;
        using Tuple3<T>::y;
        using Tuple3<T>::z;
    };

    template <typename T>
    class Tangent3 : public Vector3<T> {
    public:
        // Tangent3 Public Methods
        Tangent3() = default;
        Tangent3(T x, T y, T z) : Vector3<T>(x, y, z) {}

        template <typename U>
        explicit Tangent3(Vector3<U> v) : Vector3<T>(T(v.x), T(v.y), T(v.z)) {}
    };

    class Point3fi : public Point3<Interval> {
    public:
        //using Point3<Interval>::x;
        //using Point3<Interval>::y;
        //using Point3<Interval>::z;
        //using Point3<Interval>::operator+;
        //using Point3<Interval>::operator*;
        //using Point3<Interval>::operator*=;

        Point3fi() = default;
        Point3fi(Interval x, Interval y, Interval z) : Point3<Interval>(x, y, z) {}
        Point3fi(float x, float y, float z)
            : Point3<Interval>(Interval(x), Interval(y), Interval(z)) {}
        Point3fi(const Point3f& p)
            : Point3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        Point3fi(Point3<Interval> p) : Point3<Interval>(p) {}
        Point3fi(Point3f p, Tangent3f e)
            : Point3<Interval>(Interval::FromValueAndError(p.x, e.x),
                Interval::FromValueAndError(p.y, e.y),
                Interval::FromValueAndError(p.z, e.z)) {}

        Tangent3f Error() const { return { x.Width() / 2, y.Width() / 2, z.Width() / 2 }; }

        bool IsExact() const { return x.Width() == 0 && y.Width() == 0 && z.Width() == 0; }

        // Meh--can't seem to get these from Point3 via using declarations...
        template <typename U>
        Point3fi operator+(Tangent3<U> v) const {
            return { x + v.x, y + v.y, z + v.z };
        }
        template <typename U>
        Point3fi& operator+=(Tangent3<U> v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        Point3fi operator-() const { return { -x, -y, -z }; }
        template <typename U>
        Point3fi operator-(Point3<U> p) const {
            return { x - p.x, y - p.y, z - p.z };
        }
        template <typename U>
        Point3fi operator-(Tangent3<U> v) const {
            return { x - v.x, y - v.y, z - v.z };
        }
        template <typename U>
        Point3fi& operator-=(Tangent3<U> v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
    };

    class Tangent3fi : public Tangent3<Interval> {
    public:
        // Tangent3fi Public Methods
        //using Tangent3<Interval>::x;
        //using Tangent3<Interval>::y;
        //using Tangent3<Interval>::z;
        //using Tangent3<Interval>::operator+;
        //using Tangent3<Interval>::operator+=;
        //using Tangent3<Interval>::operator*;
        //using Tangent3<Interval>::operator*=;

        Tangent3fi() = default;
        Tangent3fi(float x, float y, float z)
            : Tangent3<Interval>(Interval(x), Interval(y), Interval(z)) {}
        Tangent3fi(Interval x, Interval y, Interval z) : Tangent3<Interval>(x, y, z) {}
        Tangent3fi(Tangent3f p)
            : Tangent3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        template <typename T>
        explicit Tangent3fi(Point3<T> p)
            : Tangent3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        Tangent3fi(Tangent3<Interval> pfi) : Tangent3<Interval>(pfi) {}
        Tangent3fi(Tangent3f v, Tangent3f e)
            : Tangent3<Interval>(Interval::FromValueAndError(v.x, e.x),
                Interval::FromValueAndError(v.y, e.y),
                Interval::FromValueAndError(v.z, e.z)) {}

        Tangent3f Error() const { return { x.Width() / 2, y.Width() / 2, z.Width() / 2 }; }

        bool IsExact() const { return x.Width() == 0 && y.Width() == 0 && z.Width() == 0; }
    };

    template <typename T>
    class Normal3 : public Vector3<T> {
    public:
        // Normal3 Public Methods
        Normal3() = default;
        Normal3(T x, T y, T z) : Vector3<T>(x, y, z) {}

        template <typename U>
        explicit Normal3<T>(Vector3<U> v) : Vector3<T>(T(v.x), T(v.y), T(v.z)) {}
    };

    // Tuple2 Inline Functions
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Abs(C<T> t) {
        return { std::abs(t.x), std::abs(t.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Ceil(C<T> t) {
        return { std::ceil(t.x), std::ceil(t.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Floor(C<T> t) {
        return { std::floor(t.x), std::floor(t.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Min(C<T> t1, C<T> t2) {
        return { std::min(t1.x, t2.x), std::min(t1.y, t2.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline T MinComponentValue(C<T> t) {
        return std::min({ t.x, t.y });
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline int MinComponentIndex(C<T> t) {
        return (t.x < t.y) ? 0 : 1;
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Max(C<T> t1, C<T> t2) {
        return { std::max(t1.x, t2.x), std::max(t1.y, t2.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline T MaxComponentValue(C<T> t) {
        return std::max({ t.x, t.y });
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline int MaxComponentIndex(C<T> t) {
        return (t.x > t.y) ? 0 : 1;
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline C<T> Permute(C<T> t, std::array<int, 2> p) {
        return { t[p[0]], t[p[1]] };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple2<T>, C<T>>::value>* = nullptr>
    inline T HProd(C<T> t) {
        return t.x * t.y;
    }

    // Tuple3 Inline Functions
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Abs(C<T> t) {
        return { std::abs(t.x), std::abs(t.y), std::abs(t.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Ceil(C<T> t) {
        return { std::ceil(t.x), std::ceil(t.y), std::ceil(t.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Floor(C<T> t) {
        return { std::floor(t.x), std::floor(t.y), std::floor(t.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Min(C<T> t1, C<T> t2) {
        return { std::min(t1.x, t2.x), std::min(t1.y, t2.y), std::min(t1.z, t2.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline T MinComponentValue(C<T> t) {
        return std::min({ t.x, t.y, t.z });
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline int MinComponentIndex(C<T> t) {
        return (t.x < t.y) ? ((t.x < t.z) ? 0 : 2) : ((t.y < t.z) ? 1 : 2);
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Max(C<T> t1, C<T> t2) {
        return { std::max(t1.x, t2.x), std::max(t1.y, t2.y), std::max(t1.z, t2.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline T MaxComponentValue(C<T> t) {
        return std::max({ t.x, t.y, t.z });
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline int MaxComponentIndex(C<T> t) {
        return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2) : ((t.y > t.z) ? 1 : 2);
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline C<T> Permute(C<T> t, std::array<int, 3> p) {
        return { t[p[0]], t[p[1]], t[p[2]] };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Tuple3<T>, C<T>>::value>* = nullptr>
    inline T HProd(C<T> t) {
        return t.x * t.y * t.z;
    }

    // Vector2 Inline Functions
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto operator+(C<T> left, C<U> right)->C<decltype(T{} + U{}) >
    {
        return { left.x + right.x,left.y + right.y };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto operator-(C<T> left, C<U> right)->C<decltype(T{} - U{}) >
    {
        return { left.x - right.x,left.y - right.y };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto operator*(C<T> left, U right)->C<decltype(T{} *U{}) > {
        return { left.x * right, left.y * right };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto operator*(U left, C<T> right)->C<decltype(T{} *U{}) > {
        return right * left;
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto operator/(C<T> left, U right)->C<decltype(T{} / U{}) > {
        return { left.x / right, left.y / right };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline C<T> operator-(C<T> t) {
        return { -t.x, -t.y };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline auto Lerp(float t, C<T> t0, C<T> t1) {
        return (1 - t) * t0 + t * t1;
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline C<T> FMA(float a, C<T> b, C<T> c) {
        return { std::fma(a, b.x, c.x), std::fma(a, b.y, c.y) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector2<T>, C<T>>::value>* = nullptr>
    inline C<T> FMA(C<T> a, float b, C<T> c) {
        return FMA(b, a, c);
    }
    template <typename T>
    inline auto Dot(Vector2<T> v1, Vector2<T> v2) -> typename TupleLength<T>::type {
        return SumOfProducts(v1.x, v2.x, v1.y, v2.y);
    }
    template <typename T>
    inline auto AbsDot(Vector2<T> v1, Vector2<T> v2) -> typename TupleLength<T>::type {
        return std::abs(Dot(v1, v2));
    }
    template <typename T>
    inline auto LengthSquared(Vector2<T> v) -> typename TupleLength<T>::type {
        return v.x * v.x + v.y * v.y;
    }
    template <typename T>
    inline auto Length(Vector2<T> v) -> typename TupleLength<T>::type {
        return std::sqrt(LengthSquared(v));
    }
    template <typename T>
    inline auto Normalize(Vector2<T> v) {
        return v / Length(v);
    }

    // Vector3 Inline Functions
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto operator+(C<T> left, C<U> right)->C<decltype(T{} + U{}) >
    {
        return { left.x + right.x,left.y + right.y,left.z + right.z };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto operator-(C<T> left, C<U> right)->C<decltype(T{} - U{}) >
    {
        return { left.x - right.x,left.y - right.y,left.z - right.z };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto operator*(C<T> left, U right)->C<decltype(T{} *U{}) > {
        return { left.x * right, left.y * right, left.z * right };
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto operator*(U left, C<T> right)->C<decltype(T{} *U{}) > {
        return right * left;
    }
    template <typename T, typename U, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto operator/(C<T> left, U right)->C<decltype(T{} / U{}) > {
        return { left.x / right, left.y / right, left.z / right };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline C<T> operator-(C<T> t) {
        return { -t.x, -t.y, -t.z };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto Lerp(float t, C<T> t0, C<T> t1) {
        return (1 - t) * t0 + t * t1;
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline C<T> FMA(float a, C<T> b, C<T> c) {
        return { std::fma(a, b.x, c.x), std::fma(a, b.y, c.y), std::fma(a, b.z, c.z) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline C<T> FMA(C<T> a, float b, C<T> c) {
        return FMA(b, a, c);
    }
    template <typename T>
    inline auto Dot(Vector3<T> v1, Vector3<T> v2) -> typename TupleLength<T>::type {
        return std::fma(v1.x, v2.x, SumOfProducts(v1.y, v2.y, v1.z, v2.z));
    }
    template <typename T>
    inline auto AbsDot(Vector3<T> v1, Vector3<T> v2) -> typename TupleLength<T>::type {
        return std::abs(Dot(v1, v2));
    }
    template <typename T>
    inline auto LengthSquared(Vector3<T> n) -> typename TupleLength<T>::type {
        return n.x * n.x + n.y * n.y + n.z * n.z;
    }
    template <typename T>
    inline auto Length(Vector3<T> v) -> typename TupleLength<T>::type {
        return std::sqrt(LengthSquared(v));
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline auto Normalize(C<T> v) {
        return v / Length(v);
    }
    template <typename T>
    inline Tangent3<T> Cross(Vector3<T> v1, Vector3<T> v2) {
        return { DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
                DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
                DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x) };
    }
    template <typename T, template <typename> class C,
        std::enable_if_t<std::is_base_of<Vector3<T>, C<T>>::value>* = nullptr>
    inline C<T> FaceForward(C<T> v1, C<T> v2) {
        return (Dot(v1, v2) < 0.f) ? -v1 : v1;
    }
    // Equivalent to std::acos(Dot(a, b)), but more numerically stable.
    // via http://www.plunk.org/~hatch/rightway.html
    template <typename T>
    inline float AngleBetween(Vector3<T> v1, Vector3<T> v2) {
        if (Dot(v1, v2) < 0)
            return Pi - 2 * SafeASin(Length(v1 + v2) / 2);
        else
            return 2 * SafeASin(Length(v2 - v1) / 2);
    }

    template <typename T>
    inline Tangent3<T> GramSchmidt(Tangent3<T> v, Tangent3<T> w) {
        return v - Dot(v, w) * w;
    }
    template <typename T>
    inline void CoordinateSystem(Tangent3<T> v1, Tangent3<T>* v2, Tangent3<T>* v3) {
        float sign = std::copysign(float(1), v1.z);
        float a = -1 / (sign + v1.z);
        float b = v1.x * v1.y * a;
        *v2 = Tangent3<T>(1 + sign * v1.x * v1.x * a, sign * b, -sign * v1.x);
        *v3 = Tangent3<T>(b, sign + v1.y * v1.y * a, -v1.y);
    }

    // Point Inline Functions
    template <typename T>
    inline Point2<T> Lerp(float t, Point2<T> t0, Point2<T> t1) {
        return { (1 - t) * t0.x + t * t1.x, (1 - t) * t0.y + t * t1.y };
    }
    template <typename T>
    inline auto Distance(Point2<T> p1, Point2<T> p2) -> typename TupleLength<T>::type {
        return Length(p1 - p2);
    }
    template <typename T>
    inline auto DistanceSquared(Point2<T> p1, Point2<T> p2) -> typename TupleLength<T>::type {
        return LengthSquared(p1 - p2);
    }
    template <typename T>
    inline Point3<T> Lerp(float t, Point3<T> t0, Point3<T> t1) {
        return { (1 - t) * t0.x + t * t1.x, (1 - t) * t0.y + t * t1.y, (1 - t) * t0.z + t * t1.z };
    }
    template <typename T>
    inline auto Distance(Point3<T> p1, Point3<T> p2) -> typename TupleLength<T>::type {
        return Length(p1 - p2);
    }
    template <typename T>
    inline auto DistanceSquared(Point3<T> p1, Point3<T> p2) -> typename TupleLength<T>::type {
        return LengthSquared(p1 - p2);
    }

    template <typename T>
    class Bounds2 {
    public:
        // Bounds2 Public Methods
        Bounds2() {
            T minNum = std::numeric_limits<T>::lowest();
            T maxNum = std::numeric_limits<T>::max();
            pMin = Point2<T>(maxNum, maxNum);
            pMax = Point2<T>(minNum, minNum);
        }
        explicit Bounds2(Point2<T> p) : pMin(p), pMax(p) {}
        Bounds2(Point2<T> p1, Point2<T> p2) : pMin(Min(p1, p2)), pMax(Max(p1, p2)) {}
        template <typename U>
        explicit Bounds2(const Bounds2<U>& b) {
            if (b.IsEmpty())
                // Be careful about overflowing float->int conversions and the
                // like.
                *this = Bounds2<T>();
            else {
                pMin = Point2<T>(b.pMin);
                pMax = Point2<T>(b.pMax);
            }
        }

        bool IsEmpty() const {
            return pMin.x >= pMax.x || pMin.y >= pMax.y;
        }
        bool IsDegenerate() const {
            return pMin.x > pMax.x || pMin.y > pMax.y;
        }

        int MaxDimension() const {
            Vector2<T> diag = Diagonal();
            if (diag.x > diag.y)
                return 0;
            else
                return 1;
        }

        Point2<T> operator[](int i) const {
            return (i == 0) ? pMin : pMax;
        }
        Point2<T>& operator[](int i) {
            return (i == 0) ? pMin : pMax;
        }

        bool operator==(const Bounds2<T>& b) const {
            return b.pMin == pMin && b.pMax == pMax;
        }
        bool operator!=(const Bounds2<T>& b) const {
            return b.pMin != pMin || b.pMax != pMax;
        }

        Vector2<T> Diagonal() const {
            return pMax - pMin;
        }
        T Area() const {
            Vector2<T> d = pMax - pMin;
            return d.x * d.y;
        }
        Point2<T> Corner(int corner) const {
            return Point2<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y);
        }
        Point2<T> Lerps(Vector2<float> t) const {
            return Point2<T>(std::lerp(t.x, pMin.x, pMax.x),
                std::lerp(t.y, pMin.y, pMax.y));
        }
        Vector2<T> Offset(Point2<T> p) const {
            Vector2<T> o = p - pMin;
            if (pMax.x > pMin.x)
                o.x /= pMax.x - pMin.x;
            if (pMax.y > pMin.y)
                o.y /= pMax.y - pMin.y;
            return o;
        }
        void BoundingSphere(Point2<T>* center, float* radius) const {
            *center = Lerp(0.5f, pMin, pMax);
            *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
        }

        // Bounds2 Public Members
        Point2<T> pMin, pMax;
    };

    template <typename T>
    class Bounds3 {
    public:
        // Bounds3 Public Methods
        Bounds3() {
            T minNum = std::numeric_limits<T>::lowest();
            T maxNum = std::numeric_limits<T>::max();
            pMin = Point3<T>(maxNum, maxNum, maxNum);
            pMax = Point3<T>(minNum, minNum, minNum);
        }
        explicit Bounds3(Point3<T> p) : pMin(p), pMax(p) {}
        Bounds3(Point3<T> p1, Point3<T> p2) : pMin(Min(p1, p2)), pMax(Max(p1, p2)) {}
        template <typename U>
        explicit Bounds3(const Bounds3<U>& b) {
            if (b.IsEmpty())
                // Be careful about overflowing float->int conversions and the
                // like.
                *this = Bounds3<T>();
            else {
                pMin = Point3<T>(b.pMin);
                pMax = Point3<T>(b.pMax);
            }
        }

        bool IsEmpty() const {
            return pMin.x >= pMax.x || pMin.y >= pMax.y || pMin.z >= pMax.z;
        }
        bool IsDegenerate() const {
            return pMin.x > pMax.x || pMin.y > pMax.y || pMin.z > pMax.z;
        }
        int MaxDimension() const {
            Tangent3<T> d = Diagonal();
            if (d.x > d.y && d.x > d.z)
                return 0;
            else if (d.y > d.z)
                return 1;
            else
                return 2;
        }

        Point3<T> operator[](int i) const {
            return (i == 0) ? pMin : pMax;
        }
        Point3<T>& operator[](int i) {
            return (i == 0) ? pMin : pMax;
        }

        bool operator==(const Bounds3<T>& b) const {
            return b.pMin == pMin && b.pMax == pMax;
        }
        bool operator!=(const Bounds3<T>& b) const {
            return b.pMin != pMin || b.pMax != pMax;
        }

        Tangent3<T> Diagonal() const {
            return pMax - pMin;
        }
        T SurfaceArea() const {
            Tangent3<T> d = Diagonal();
            return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
        }
        T Volume() const {
            Tangent3<T> d = Diagonal();
            return d.x * d.y * d.z;
        }
        Point3<T> Corner(int corner) const {
            return Point3<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y,
                (*this)[(corner & 4) ? 1 : 0].z);
        }
        Point3f Lerps(Vector3<float> t) const {
            return Point3f(std::lerp(t.x, pMin.x, pMax.x), std::lerp(t.y, pMin.y, pMax.y),
                std::lerp(t.z, pMin.z, pMax.z));
        }
        Tangent3f Offset(Point3f p) const {
            Tangent3f o = p - pMin;
            if (pMax.x > pMin.x)
                o.x /= pMax.x - pMin.x;
            if (pMax.y > pMin.y)
                o.y /= pMax.y - pMin.y;
            if (pMax.z > pMin.z)
                o.z /= pMax.z - pMin.z;
            return o;
        }
        void BoundingSphere(Point3<T>* center, float* radius) const {
            *center = Lerp(0.5f, pMin, pMax);
            *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
        }

        bool IntersectP(Point3f o, Tangent3f d, float tMax = Infinity, float* hitt0 = nullptr,
            float* hitt1 = nullptr) const;
        bool IntersectP(Point3f o, Tangent3f d, float tMax, Tangent3f invDir,
            const int dirIsNeg[3]) const;

        // Bounds3 Public Members
        Point3<T> pMin, pMax;
    };

    class Bounds2iIterator : public std::forward_iterator_tag {
    public:
        Bounds2iIterator(const Bounds2i& b, const Point2i& pt) : p(pt), bounds(&b) {}
        Bounds2iIterator operator++() {
            advance();
            return *this;
        }
        Bounds2iIterator operator++(int) {
            Bounds2iIterator old = *this;
            advance();
            return old;
        }

        bool operator==(const Bounds2iIterator& bi) const {
            return p == bi.p && bounds == bi.bounds;
        }
        bool operator!=(const Bounds2iIterator& bi) const {
            return p != bi.p || bounds != bi.bounds;
        }

        Point2i operator*() const {
            return p;
        }

    private:
        void advance() {
            ++p.x;
            if (p.x == bounds->pMax.x) {
                p.x = bounds->pMin.x;
                ++p.y;
            }
        }
        Point2i p;
        const Bounds2i* bounds;
    };

    // Bounds2 Inline Functions
    template <typename T>
    inline Bounds2<T> Union(const Bounds2<T>& b, Point2<T> p) {
        // Be careful to not run the two-point Bounds constructor.
        Bounds2<T> ret;
        ret.pMin = Min(b.pMin, p);
        ret.pMax = Max(b.pMax, p);
        return ret;
    }
    template <typename T>
    inline Bounds2<T> Union(const Bounds2<T>& b1, const Bounds2<T>& b2) {
        // Be careful to not run the two-point Bounds constructor.
        Bounds2<T> ret;
        ret.pMin = Min(b1.pMin, b2.pMin);
        ret.pMax = Max(b1.pMax, b2.pMax);
        return ret;
    }
    template <typename T>
    inline Bounds2<T> Intersect(const Bounds2<T>& b1, const Bounds2<T>& b2) {
        // Be careful to not run the two-point Bounds constructor.
        Bounds2<T> b;
        b.pMin = Max(b1.pMin, b2.pMin);
        b.pMax = Min(b1.pMax, b2.pMax);
        return b;
    }
    template <typename T>
    inline bool Overlaps(const Bounds2<T>& ba, const Bounds2<T>& bb) {
        bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
        bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
        return (x && y);
    }
    template <typename T>
    inline bool Inside(Point2<T> pt, const Bounds2<T>& b) {
        return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y && pt.y <= b.pMax.y);
    }
    template <typename T>
    inline bool Inside(const Bounds2<T>& ba, const Bounds2<T>& bb) {
        return (ba.pMin.x >= bb.pMin.x && ba.pMax.x <= bb.pMax.x && ba.pMin.y >= bb.pMin.y &&
            ba.pMax.y <= bb.pMax.y);
    }
    template <typename T>
    inline bool InsideExclusive(Point2<T> pt, const Bounds2<T>& b) {
        return (pt.x >= b.pMin.x && pt.x < b.pMax.x&& pt.y >= b.pMin.y && pt.y < b.pMax.y);
    }
    template <typename T, typename U>
    inline Bounds2<T> Expand(const Bounds2<T>& b, U delta) {
        Bounds2<T> ret;
        ret.pMin = b.pMin - Vector2<T>(delta, delta);
        ret.pMax = b.pMax + Vector2<T>(delta, delta);
        return ret;
    }

    // Bounds3 Inline Functions
    template <typename T>
    inline Bounds3<T> Union(const Bounds3<T>& b, Point3<T> p) {
        Bounds3<T> ret;
        ret.pMin = Min(b.pMin, p);
        ret.pMax = Max(b.pMax, p);
        return ret;
    }
    template <typename T>
    inline Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2) {
        Bounds3<T> ret;
        ret.pMin = Min(b1.pMin, b2.pMin);
        ret.pMax = Max(b1.pMax, b2.pMax);
        return ret;
    }
    template <typename T>
    inline Bounds3<T> Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2) {
        Bounds3<T> b;
        b.pMin = Max(b1.pMin, b2.pMin);
        b.pMax = Min(b1.pMax, b2.pMax);
        return b;
    }
    template <typename T>
    inline bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2) {
        bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
        bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
        bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
        return (x && y && z);
    }
    template <typename T>
    inline bool Inside(Point3<T> p, const Bounds3<T>& b) {
        return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y && p.y <= b.pMax.y &&
            p.z >= b.pMin.z && p.z <= b.pMax.z);
    }
    template <typename T>
    inline bool InsideExclusive(Point3<T> p, const Bounds3<T>& b) {
        return (p.x >= b.pMin.x && p.x < b.pMax.x&& p.y >= b.pMin.y && p.y < b.pMax.y&&
            p.z >= b.pMin.z && p.z < b.pMax.z);
    }
    template <typename T, typename U>
    inline auto DistanceSquared(Point3<T> p, const Bounds3<U>& b) {
        using TDist = decltype(T{} - U{});
        TDist dx = std::max<TDist>({ 0, b.pMin.x - p.x, p.x - b.pMax.x });
        TDist dy = std::max<TDist>({ 0, b.pMin.y - p.y, p.y - b.pMax.y });
        TDist dz = std::max<TDist>({ 0, b.pMin.z - p.z, p.z - b.pMax.z });
        return Sqr(dx) + Sqr(dy) + Sqr(dz);
    }
    template <typename T, typename U>
    inline auto Distance(Point3<T> p, const Bounds3<U>& b) {
        auto dist2 = DistanceSquared(p, b);
        using TDist = typename TupleLength<decltype(dist2)>::type;
        return std::sqrt(TDist(dist2));
    }
    template <typename T, typename U>
    inline Bounds3<T> Expand(const Bounds3<T>& b, U delta) {
        Bounds3<T> ret;
        ret.pMin = b.pMin - Tangent3<T>(delta, delta, delta);
        ret.pMax = b.pMax + Tangent3<T>(delta, delta, delta);
        return ret;
    }
    template <typename T>
    inline bool Bounds3<T>::IntersectP(Point3f o, Tangent3f d, float tMax,
        float* hitt0, float* hitt1) const {
        float t0 = 0, t1 = tMax;
        for (int i = 0; i < 3; ++i) {
            // Update interval for _i_th bounding box slab
            float invRayDir = 1 / d[i];
            float tNear = (pMin[i] - o[i]) * invRayDir;
            float tFar = (pMax[i] - o[i]) * invRayDir;
            // Update parametric interval from slab intersection $t$ values
            if (tNear > tFar)
                std::swap(tNear, tFar);
            // Update _tFar_ to ensure robust ray--bounds intersection
            tFar *= 1 + 2 * gamma(3);

            t0 = tNear > t0 ? tNear : t0;
            t1 = tFar < t1 ? tFar : t1;
            if (t0 > t1)
                return false;
        }
        if (hitt0)
            *hitt0 = t0;
        if (hitt1)
            *hitt1 = t1;
        return true;
    }
    template <typename T>
    inline bool Bounds3<T>::IntersectP(Point3f o, Tangent3f d, float raytMax,
        Tangent3f invDir,
        const int dirIsNeg[3]) const {
        const Bounds3f& bounds = *this;
        // Check for ray intersection against $x$ and $y$ slabs
        float tMin = (bounds[dirIsNeg[0]].x - o.x) * invDir.x;
        float tMax = (bounds[1 - dirIsNeg[0]].x - o.x) * invDir.x;
        float tyMin = (bounds[dirIsNeg[1]].y - o.y) * invDir.y;
        float tyMax = (bounds[1 - dirIsNeg[1]].y - o.y) * invDir.y;
        // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
        tMax *= 1 + 2 * gamma(3);
        tyMax *= 1 + 2 * gamma(3);

        if (tMin > tyMax || tyMin > tMax)
            return false;
        if (tyMin > tMin)
            tMin = tyMin;
        if (tyMax < tMax)
            tMax = tyMax;

        // Check for ray intersection against $z$ slab
        float tzMin = (bounds[dirIsNeg[2]].z - o.z) * invDir.z;
        float tzMax = (bounds[1 - dirIsNeg[2]].z - o.z) * invDir.z;
        // Update _tzMax_ to ensure robust bounds intersection
        tzMax *= 1 + 2 * gamma(3);

        if (tMin > tzMax || tzMin > tMax)
            return false;
        if (tzMin > tMin)
            tMin = tzMin;
        if (tzMax < tMax)
            tMax = tzMax;

        return (tMin < raytMax) && (tMax > 0);
    }

    inline Bounds2iIterator begin(const Bounds2i& b) {
        return Bounds2iIterator(b, b.pMin);
    }
    inline Bounds2iIterator end(const Bounds2i& b) {
        // Normally, the ending point is at the minimum x value and one past
        // the last valid y value.
        Point2i pEnd(b.pMin.x, b.pMax.y);
        // However, if the bounds are degenerate, override the end point to
        // equal the start point so that any attempt to iterate over the bounds
        // exits out immediately.
        if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
            pEnd = b.pMin;
        return Bounds2iIterator(b, pEnd);
    }

    class Frame {
    public:
        // Frame Public Methods
        Frame() : x(1, 0, 0), y(0, 1, 0), z(0, 0, 1) {}
        Frame(Tangent3f x, Tangent3f y, Tangent3f z) : x(x), y(y), z(z) {}

        static Frame FromXZ(Tangent3f x, Tangent3f z) {
            return Frame(x, Cross(z, x), z);
        }
        static Frame FromXY(Tangent3f x, Tangent3f y) {
            return Frame(x, y, Cross(x, y));
        }
        static Frame FromX(Tangent3f x) {
            Tangent3f y, z;
            CoordinateSystem(x, &y, &z);
            return Frame(x, y, z);
        }
        static Frame FromY(Tangent3f y) {
            Tangent3f x, z;
            CoordinateSystem(y, &z, &x);
            return Frame(x, y, z);
        }
        static Frame FromZ(Tangent3f z) {
            Tangent3f x, y;
            CoordinateSystem(z, &x, &y);
            return Frame(x, y, z);
        }
        static Frame FromX(Normal3f x) {
            return FromX((Tangent3f)x);
        }
        static Frame FromY(Normal3f y) {
            return FromY((Tangent3f)y);
        }
        static Frame FromZ(Normal3f z) {
            return FromZ((Tangent3f)z);
        }

        Tangent3f ToLocal(Tangent3f v) const {
            return Tangent3f(Dot(v, x), Dot(v, y), Dot(v, z));
        }
        Normal3f ToLocal(Normal3f n) const {
            return Normal3f(Dot(n, x), Dot(n, y), Dot(n, z));
        }
        Tangent3f FromLocal(Tangent3f v) const {
            return v.x * x + v.y * y + v.z * z;
        }
        Normal3f FromLocal(Normal3f n) const {
            return (Normal3f)(n.x * x + n.y * y + n.z * z);
        }

        // Frame Public Members
        Tangent3f x, y, z;
    };

    // Spherical Geometry Inline Functions
    inline float SphericalTriangleArea(Tangent3f a, Tangent3f b, Tangent3f c) {
        return std::abs(
            2 * std::atan2(Dot(a, Cross(b, c)), 1 + Dot(a, b) + Dot(a, c) + Dot(b, c)));
    }
    inline float SphericalQuadArea(Tangent3f a, Tangent3f b, Tangent3f c, Tangent3f d) {
        Tangent3f axb = Cross(a, b), bxc = Cross(b, c),
            cxd = Cross(c, d), dxa = Cross(d, a);
        if (LengthSquared(axb) == 0 || LengthSquared(bxc) == 0 || LengthSquared(cxd) == 0 ||
            LengthSquared(dxa) == 0)
            return 0;
        axb = Normalize(axb);
        bxc = Normalize(bxc);
        cxd = Normalize(cxd);
        dxa = Normalize(dxa);

        float alpha = AngleBetween(dxa, -axb);
        float beta = AngleBetween(axb, -bxc);
        float gamma = AngleBetween(bxc, -cxd);
        float delta = AngleBetween(cxd, -dxa);

        return std::abs(alpha + beta + gamma + delta - 2 * Pi);
    }
    inline Vector3f SphericalDirection(float sinTheta, float cosTheta, float phi) {
        return Vector3f(Clamp(sinTheta, -1, 1) * std::cos(phi),
            Clamp(sinTheta, -1, 1) * std::sin(phi), Clamp(cosTheta, -1, 1));
    }
    inline float SphericalTheta(Vector3f v) {
        return SafeACos(v.z);
    }
    inline float SphericalPhi(Vector3f v) {
        float p = std::atan2(v.y, v.x);
        return (p < 0) ? (p + 2 * Pi) : p;
    }
    inline float CosTheta(Vector3f w) {
        return w.z;
    }
    inline float Cos2Theta(Vector3f w) {
        return w.z * w.z;
    }
    inline float AbsCosTheta(Vector3f w) {
        return std::abs(w.z);
    }
    inline float Sin2Theta(Vector3f w) {
        return std::max<float>(0, 1 - Cos2Theta(w));
    }
    inline float SinTheta(Vector3f w) {
        return std::sqrt(Sin2Theta(w));
    }
    inline float TanTheta(Vector3f w) {
        return SinTheta(w) / CosTheta(w);
    }
    inline float Tan2Theta(Vector3f w) {
        return Sin2Theta(w) / Cos2Theta(w);
    }
    inline float CosPhi(Vector3f w) {
        float sinTheta = SinTheta(w);
        return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
    }
    inline float SinPhi(Vector3f w) {
        float sinTheta = SinTheta(w);
        return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
    }
    inline float CosDPhi(Vector3f wa, Vector3f wb) {
        float waxy = wa.x * wa.x + wa.y * wa.y, wbxy = wb.x * wb.x + wb.y * wb.y;
        if (waxy == 0 || wbxy == 0)
            return 1;
        return Clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt(waxy * wbxy), -1, 1);
    }
    inline bool SameHemisphere(Vector3f w, Vector3f wp) {
        return w.z * wp.z > 0;
    }

    class DirectionCone {
    public:
        // DirectionCone Public Methods
        DirectionCone() = default;
        DirectionCone(Tangent3f w, float cosTheta) : w(Normalize(w)), cosTheta(cosTheta) {}
        explicit DirectionCone(Tangent3f w) : DirectionCone(w, 1) {}

        bool IsEmpty() const {
            return cosTheta == Infinity;
        }
        static DirectionCone EntireSphere() {
            return DirectionCone(Tangent3f(0, 0, 1), -1);
        }
        Tangent3f ClosestVectorInCone(Tangent3f wp) const;

        // DirectionCone Public Members
        Tangent3f w;
        float cosTheta = Infinity;
    };

    // DirectionCone Function Declarations
    DirectionCone Union(const DirectionCone& a, const DirectionCone& b);

    // DirectionCone Inline Functions
    inline bool Inside(const DirectionCone& d, Tangent3f w) {
        return !d.IsEmpty() && Dot(d.w, Normalize(w)) >= d.cosTheta;
    }

    inline DirectionCone BoundSubtendedDirections(const Bounds3f& b, Point3f p) {
        // Compute bounding sphere for _b_ and check if _p_ is inside
        float radius;
        Point3f pCenter;
        b.BoundingSphere(&pCenter, &radius);
        if (DistanceSquared(p, pCenter) < radius * radius)
            return DirectionCone::EntireSphere();
        // Compute and return _DirectionCone_ for bounding sphere
        Tangent3f w = Normalize(pCenter - p);
        float sin2ThetaMax = radius * radius / DistanceSquared(pCenter, p);
        float cosThetaMax = SafeSqrt(1 - sin2ThetaMax);
        return DirectionCone(w, cosThetaMax);
    }
    inline Tangent3f DirectionCone::ClosestVectorInCone(Tangent3f wp) const {
        wp = Normalize(wp);
        // Return provided vector if it is inside the cone
        if (Dot(wp, w) > cosTheta)
            return wp;
        // Find closest vector by rotating _wp_ until it touches the cone
        float sinTheta = -SafeSqrt(1 - cosTheta * cosTheta);
        Tangent3f a = Cross(wp, w);
        return cosTheta * w +
            sinTheta / Length(a) *
            Tangent3f(w.x * (wp.y * w.y + wp.z * w.z) - wp.x * (w.y * w.y + w.z * w.z),
                w.y * (wp.x * w.x + wp.z * w.z) - wp.y * (w.x * w.x + w.z * w.z),
                w.z * (wp.x * w.x + wp.y * w.y) - wp.z * (w.x * w.x + w.y * w.y));
    }

    inline Point3f OffsetRayOrigin(const Point3f& p, const Tangent3f& pError,
        const Normal3f& n, const Tangent3f& w) {
        float d = Dot(Abs(n), pError);
        Tangent3f offset = d * Tangent3f(n);
        if (Dot(w, n) < 0) offset = -offset;
        Point3f po = p + offset;
        // Round offset point _po_ away from _p_
        for (int i = 0; i < 3; ++i) {
            if (offset[i] > 0)
                po[i] = NextFloatUp(po[i]);
            else if (offset[i] < 0)
                po[i] = NextFloatDown(po[i]);
        }
        return po;
    }

} // namespace lightfold