#pragma once

#include <array>

#include <math/math.h>

namespace lightfold {

    template <typename T>
    class GeoVector2;
    template <typename T>
    class GeoVector3;
    using GeoVector2f = GeoVector2<float>;
    using GeoVector2i = GeoVector2<int>;
    using GeoVector3f = GeoVector3<float>;
    using GeoVector3i = GeoVector3<int>;

    template <typename T>
    class Point2;
    template <typename T>
    class Point3;
    using Point2f = Point2<float>;
    using Point2i = Point2<int>;
    using Point3f = Point3<float>;
    using Point3i = Point3<int>;

    template <typename T>
    class TanVector3;
    using TanVector3f = TanVector3<float>;

    template <typename T>
    class CotVector3;
    using CotVector3f = CotVector3<float>;

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
    }  // anonymous namespace

    // Class Definitions
    template <template <typename> class Rep, typename T>
    class Tuple2 {
    public:
        static const int nDimensions = 2;

        // Tuple2 Public Methods
        Tuple2() = default;
        Tuple2(T x, T y) : x(x), y(y) {}

        template <typename U>
        auto operator+(Rep<U> c) const -> Rep<decltype(T{} + U{})> {
            return { x + c.x, y + c.y };
        }
        template <typename U>
        Rep<T>& operator+=(Rep<U> c) {
            x += c.x;
            y += c.y;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator-(Rep<U> c) const -> Rep<decltype(T{} - U{})> {
            return { x - c.x, y - c.y };
        }
        template <typename U>
        Rep<T>& operator-=(Rep<U> c) {
            x -= c.x;
            y -= c.y;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator*(U s) const -> Rep<decltype(T{} * U{})> {
            return { s * x, s * y };
        }
        template <typename U>
        Rep<T>& operator*=(U s) {
            x *= s;
            y *= s;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator/(U d) const -> Rep<decltype(T{} / U{})> {
            return { x / d, y / d };
        }
        template <typename U>
        Rep<T>& operator/=(U d) {
            x /= d;
            y /= d;
            return static_cast<Rep<T> &>(*this);
        }
        Rep<T> operator-() const {
            return { -x, -y };
        }

        bool operator==(Rep<T> c) const {
            return x == c.x && y == c.y;
        }
        bool operator!=(Rep<T> c) const {
            return x != c.x || y != c.y;
        }

        T operator[](int i) const {
            return (i == 0) ? x : y;
        }
        T& operator[](int i) {
            return (i == 0) ? x : y;
        }

        // Tuple2 Public Members
        T x{}, y{};
    };

    template <template <typename> class Rep, typename T>
    class Tuple3 {
    public:
        static const int nDimensions = 3;

        // Tuple3 Public Methods
        Tuple3() = default;
        Tuple3(T x, T y, T z) : x(x), y(y), z(z) { }

        template <typename U>
        auto operator+(Rep<U> c) const -> Rep<decltype(T{} + U{})> {
            return { x + c.x, y + c.y, z + c.z };
        }
        template <typename U>
        Rep<T>& operator+=(Rep<U> c) {
            x += c.x;
            y += c.y;
            z += c.z;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator-(Rep<U> c) const -> Rep<decltype(T{} - U{})> {
            return { x - c.x, y - c.y, z - c.z };
        }
        template <typename U>
        Rep<T>& operator-=(Rep<U> c) {
            x -= c.x;
            y -= c.y;
            z -= c.z;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator*(U s) const -> Rep<decltype(T{} * U{})> {
            return { s * x, s * y, s * z };
        }
        template <typename U>
        Rep<T>& operator*=(U s) {
            x *= s;
            y *= s;
            z *= s;
            return static_cast<Rep<T> &>(*this);
        }
        template <typename U>
        auto operator/(U d) const -> Rep<decltype(T{} / U{}) > {
            return { x / d, y / d, z / d };
        }
        template <typename U>
        Rep<T>& operator/=(U d) {
            x /= d;
            y /= d;
            z /= d;
            return static_cast<Rep<T> &>(*this);
        }
        Rep<T> operator-() const {
            return { -x, -y, -z };
        }

        bool operator==(Rep<T> c) const {
            return x == c.x && y == c.y && z == c.z;
        }
        bool operator!=(Rep<T> c) const {
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
        T x{}, y{}, z{};
    };

    template <typename T>
    class Point2 : public Tuple2<Point2, T> {
    public:
        // Point2 Public Methods
        //using Tuple2<Point2, T>::x;
        //using Tuple2<Point2, T>::y;
        using Tuple2<Point2, T>::operator+;
        using Tuple2<Point2, T>::operator+=;
        using Tuple2<Point2, T>::operator*;
        using Tuple2<Point2, T>::operator*=;

        Point2() { x = y = 0; }
        Point2(T x, T y) : Tuple2<lightfold::Point2, T>(x, y) {}

        template <typename U>
        explicit Point2(Point2<U> v) : Tuple2<lightfold::Point2, T>(T(v.x), T(v.y)) {}
        template <typename U>
        explicit Point2(GeoVector2<U> v) : Tuple2<lightfold::Point2, T>(T(v.x), T(v.y)) {}

        template <typename U>
        auto operator+(GeoVector2<U> v) const -> Point2<decltype(T{} + U{})> {
            return { x + v.x, y + v.y };
        }
        template <typename U>
        Point2<T>& operator+=(GeoVector2<U> v) {
            x += v.x;
            y += v.y;
            return *this;
        }
        // We can't do using operator- above, since we don't want to pull in
        // the Point-Point -> Point one so that we can return a vector
        // instead...
        Point2<T> operator-() const {
            return { -x, -y };
        }
        template <typename U>
        auto operator-(Point2<U> p) const -> GeoVector2<decltype(T{} - U{})> {
            return { x - p.x, y - p.y };
        }
        template <typename U>
        auto operator-(GeoVector2<U> v) const -> Point2<decltype(T{} - U{})> {
            return { x - v.x, y - v.y };
        }
        template <typename U>
        Point2<T>& operator-=(GeoVector2<U> v) {
            x -= v.x;
            y -= v.y;
            return *this;
        }
    };

    template <typename T>
    class Point3 : public Tuple3<Point3, T> {
    public:
        // Point3 Public Methods
        //using Tuple3<Point3, T>::x;
        //using Tuple3<Point3, T>::y;
        //using Tuple3<Point3, T>::z;
        using Tuple3<Point3, T>::operator+;
        using Tuple3<Point3, T>::operator+=;
        using Tuple3<Point3, T>::operator*;
        using Tuple3<Point3, T>::operator*=;

        Point3() = default;
        Point3(T x, T y, T z) : Tuple3<lightfold::Point3, T>(x, y, z) {}

        template <typename U>
        explicit Point3(Point3<U> p) : Tuple3<lightfold::Point3, T>(T(p.x), T(p.y), T(p.z)) {}
        template <typename U>
        explicit Point3(GeoVector3<U> v) : Tuple3<lightfold::Point3, T>(T(v.x), T(v.y), T(v.z)) {}

        template <typename U>
        auto operator+(TanVector3<U> v) const -> Point3<decltype(T{} + U{})> {
            return { x + v.x, y + v.y, z + v.z };
        }
        template <typename U>
        Point3<T>& operator+=(TanVector3<U> v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        // We can't do using operator- above, since we don't want to pull in
        // the Point-Point -> Point one so that we can return a vector
        // instead...
        Point3<T> operator-() const {
            return { -x, -y, -z };
        }
        template <typename U>
        auto operator-(Point3<U> p) const -> TanVector3<decltype(T{} - U{})> {
            return { x - p.x, y - p.y, z - p.z };
        }
        template <typename U>
        auto operator-(TanVector3<U> v) const -> Point3<decltype(T{} - U{})> {
            return { x - v.x, y - v.y, z - v.z };
        }
        template <typename U>
        Point3<T>& operator-=(TanVector3<U> v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
    };

    /*class Point3fi : public Point3<Interval> {
    public:
        using Point3<Interval>::x;
        using Point3<Interval>::y;
        using Point3<Interval>::z;
        using Point3<Interval>::operator+;
        using Point3<Interval>::operator*;
        using Point3<Interval>::operator*=;

        Point3fi() = default;
        Point3fi(Interval x, Interval y, Interval z) : Point3<Interval>(x, y, z) {}
        Point3fi(float x, float y, float z) : Point3<Interval>(Interval(x), Interval(y), Interval(z)) {}
        Point3fi(const Point3f& p) : Point3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        Point3fi(Point3<Interval> p) : Point3<Interval>(p) {}
        Point3fi(Point3f p, Vector3f e)
            : Point3<Interval>(Interval::FromValueAndError(p.x, e.x),
                Interval::FromValueAndError(p.y, e.y),
                Interval::FromValueAndError(p.z, e.z)) {}
        
        Vector3f Error() const { return { x.Width() / 2, y.Width() / 2, z.Width() / 2 }; }
        bool IsExact() const { return x.Width() == 0 && y.Width() == 0 && z.Width() == 0; }

        // Meh--can't seem to get these from Point3 via using declarations...
        template <typename U>
        Point3fi operator+(Vector3<U> v) const {
            return { x + v.x, y + v.y, z + v.z };
        }
        template <typename U>
        Point3fi& operator+=(Vector3<U> v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        Point3fi operator-() const {
            return { -x, -y, -z };
        }
        template <typename U>
        Point3fi operator-(Point3<U> p) const { return { x - p.x, y - p.y, z - p.z }; }
        template <typename U>
        Point3fi operator-(Vector3<U> v) const { return { x - v.x, y - v.y, z - v.z }; }
        template <typename U>
        Point3fi& operator-=(Vector3<U> v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
    };*/

    template <typename T>
    class GeoVector2 : public Tuple2<GeoVector2, T> {
    public:
        // GeoVector2 Public Methods
        //using Tuple2<GeoVector2, T>::x;
        //using Tuple2<GeoVector2, T>::y;

        GeoVector2() = default;
        GeoVector2(T x, T y) : Tuple2<lightfold::GeoVector2, T>(x, y) {}

        template <typename U>
        explicit GeoVector2(GeoVector2<U> v) : Tuple2<lightfold::GeoVector2, T>(T(v.x), T(v.y)) {}
        template <typename U>
        explicit GeoVector2(Point2<U> p) : Tuple2<lightfold::GeoVector2, T>(T(p.x), T(p.y)) {}
    };

    template <typename T>
    class GeoVector3 : public Tuple3<GeoVector3, T> {
    public:
        // Vector3 Public Methods
        //using Tuple3<GeoVector3, T>::x;
        //using Tuple3<GeoVector3, T>::y;
        //using Tuple3<GeoVector3, T>::z;

        GeoVector3() = default;
        GeoVector3(T x, T y, T z) : Tuple3<lightfold::GeoVector3, T>(x, y, z) {}

        template <typename U>
        explicit GeoVector3(GeoVector3<U> v) : Tuple3<lightfold::GeoVector3, T>(T(v.x), T(v.y), T(v.z)) {}
        template <typename U>
        explicit GeoVector3(Point3<U> p) : Tuple3<lightfold::GeoVector3, T>(T(p.x), T(p.y), T(p.z)) {}
        template <typename U>
        explicit GeoVector3(TanVector3<U> n) : Tuple3<lightfold::GeoVector3, T>(T(n.x), T(n.y), T(n.z)) {}
        template <typename U>
        explicit GeoVector3(CotVector3<U> n) : Tuple3<lightfold::GeoVector3, T>(T(n.x), T(n.y), T(n.z)) {}
    };

    /*class Vector3fi : public Vector3<Interval> {
    public:
        // Vector3fi Public Methods
        using Vector3<Interval>::x;
        using Vector3<Interval>::y;
        using Vector3<Interval>::z;
        using Vector3<Interval>::operator+;
        using Vector3<Interval>::operator+=;
        using Vector3<Interval>::operator*;
        using Vector3<Interval>::operator*=;

        Vector3fi() = default;
        Vector3fi(float x, float y, float z)  : Vector3<Interval>(Interval(x), Interval(y), Interval(z)) {}
        Vector3fi(Interval x, Interval y, Interval z) : Vector3<Interval>(x, y, z) {}
        Vector3fi(Vector3f p) : Vector3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        template <typename T>
        explicit Vector3fi(Point3<T> p) : Vector3<Interval>(Interval(p.x), Interval(p.y), Interval(p.z)) {}
        Vector3fi(Vector3<Interval> pfi) : Vector3<Interval>(pfi) {}
        Vector3fi(Vector3f v, Vector3f e)
            : Vector3<Interval>(Interval::FromValueAndError(v.x, e.x),
                Interval::FromValueAndError(v.y, e.y),
                Interval::FromValueAndError(v.z, e.z)) {}
        Vector3f Error() const { return { x.Width() / 2, y.Width() / 2, z.Width() / 2 }; }
        bool IsExact() const { return x.Width() == 0 && y.Width() == 0 && z.Width() == 0; }
    };*/

    template <typename T>
    class TanVector3 : public GeoVector3<T> {
    public:
        // TanVector3 Public Methods
        using GeoVector3<T>::operator+;
        using GeoVector3<T>::operator+=;
        using GeoVector3<T>::operator-;
        using GeoVector3<T>::operator-=;
        using GeoVector3<T>::operator*;
        using GeoVector3<T>::operator*=;
        using GeoVector3<T>::operator/;
        using GeoVector3<T>::operator/=;
        using GeoVector3<T>::operator==;
        using GeoVector3<T>::operator!=;
        using GeoVector3<T>::operator[];

        TanVector3() = default;
        TanVector3(T x, T y, T z) : GeoVector3<T>(x, y, z) {}

        template <typename U>
        explicit TanVector3(GeoVector3<U> v) : GeoVector3<T>(T(v.x), T(v.y), T(v.z)) {}
    };

    template <typename T>
    class CotVector3 : public GeoVector3<T> {
    public:
        // CotVector3 Public Methods
        using GeoVector3<T>::operator+;
        using GeoVector3<T>::operator+=;
        using GeoVector3<T>::operator-;
        using GeoVector3<T>::operator-=;
        using GeoVector3<T>::operator*;
        using GeoVector3<T>::operator*=;
        using GeoVector3<T>::operator/;
        using GeoVector3<T>::operator/=;
        using GeoVector3<T>::operator==;
        using GeoVector3<T>::operator!=;
        using GeoVector3<T>::operator[];

        CotVector3() = default;
        CotVector3(T x, T y, T z) : GeoVector3<T>(x, y, z) {}

        template <typename U>
        explicit CotVector3<T>(GeoVector3<U> v) : GeoVector3<T>(T(v.x), T(v.y), T(v.z)) {}
    };

    class Quaternion {
    public:
        // Quaternion Public Methods
        Quaternion() = default;

        Quaternion& operator+=(Quaternion q) {
            v += q.v;
            w += q.w;
            return *this;
        }
        Quaternion operator+(Quaternion q) const {
            return { v + q.v, w + q.w };
        }
        Quaternion& operator-=(Quaternion q) {
            v -= q.v;
            w -= q.w;
            return *this;
        }
        Quaternion operator-() const {
            return { -v, -w };
        }
        Quaternion operator-(Quaternion q) const {
            return { v - q.v, w - q.w };
        }
        Quaternion& operator*=(float f) {
            v *= f;
            w *= f;
            return *this;
        }
        Quaternion operator*(float f) const {
            return { v * f, w * f };
        }
        Quaternion& operator/=(float f) {
            v /= f;
            w /= f;
            return *this;
        }
        Quaternion operator/(float f) const {
            return { v / f, w / f };
        }

        // Quaternion Public Members
        GeoVector3f v;
        float w = 1;
    };

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

        GeoVector2<T> Diagonal() const {
            return pMax - pMin;
        }
        T Area() const {
            GeoVector2<T> d = pMax - pMin;
            return d.x * d.y;
        }

        bool IsEmpty() const {
            return pMin.x >= pMax.x || pMin.y >= pMax.y;
        }
        bool IsDegenerate() const {
            return pMin.x > pMax.x || pMin.y > pMax.y;
        }

        int MaxDimension() const {
            GeoVector2<T> diag = Diagonal();
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

        Point2<T> Corner(int corner) const {
            return Point2<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y);
        }
        Point2<T> Lerp(Point2f t) const {
            return Point2<T>(lightfold::Lerp(t.x, pMin.x, pMax.x),
                lightfold::Lerp(t.y, pMin.y, pMax.y));
        }
        GeoVector2<T> Offset(Point2<T> p) const {
            GeoVector2<T> o = p - pMin;
            if (pMax.x > pMin.x)
                o.x /= pMax.x - pMin.x;
            if (pMax.y > pMin.y)
                o.y /= pMax.y - pMin.y;
            return o;
        }
        void BoundingSphere(Point2<T>* c, float* rad) const {
            *c = (pMin + pMax) / 2;
            *rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
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

        Point3<T> operator[](int i) const {
            return (i == 0) ? pMin : pMax;
        }
        Point3<T>& operator[](int i) {
            return (i == 0) ? pMin : pMax;
        }

        Point3<T> Corner(int corner) const {
            return Point3<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y,
                (*this)[(corner & 4) ? 1 : 0].z);
        }
        TanVector3<T> Diagonal() const {
            return pMax - pMin;
        }
        T SurfaceArea() const {
            Vector3<T> d = Diagonal();
            return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
        }
        T Volume() const {
            Vector3<T> d = Diagonal();
            return d.x * d.y * d.z;
        }
        int MaxDimension() const {
            Vector3<T> d = Diagonal();
            if (d.x > d.y && d.x > d.z)
                return 0;
            else if (d.y > d.z)
                return 1;
            else
                return 2;
        }

        Point3f Lerp(Point3f t) const {
            return Point3f(lightfold::Lerp(t.x, pMin.x, pMax.x), lightfold::Lerp(t.y, pMin.y, pMax.y),
                lightfold::Lerp(t.z, pMin.z, pMax.z));
        }
        TanVector3f Offset(Point3f p) const {
            TanVector3f o = p - pMin;
            if (pMax.x > pMin.x)
                o.x /= pMax.x - pMin.x;
            if (pMax.y > pMin.y)
                o.y /= pMax.y - pMin.y;
            if (pMax.z > pMin.z)
                o.z /= pMax.z - pMin.z;
            return o;
        }
        void BoundingSphere(Point3<T>* center, float* radius) const {
            *center = (pMin + pMax) / 2;
            *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
        }
        bool IsEmpty() const {
            return pMin.x >= pMax.x || pMin.y >= pMax.y || pMin.z >= pMax.z;
        }
        bool IsDegenerate() const {
            return pMin.x > pMax.x || pMin.y > pMax.y || pMin.z > pMax.z;
        }

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
        bool operator==(const Bounds3<T>& b) const {
            return b.pMin == pMin && b.pMax == pMax;
        }
        bool operator!=(const Bounds3<T>& b) const {
            return b.pMin != pMin || b.pMax != pMax;
        }

        bool IntersectP(Point3f o, TanVector3f d, float tMax = Infinity, float* hitt0 = nullptr,
            float* hitt1 = nullptr) const;
        bool IntersectP(Point3f o, TanVector3f d, float tMax, TanVector3f invDir,
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

    // Tuple3 Inline Functions
    template <template <class> class C, typename T, typename U>
    inline auto operator*(U s, Tuple3<C, T> t)->C<decltype(T{} *U{}) > {
        return t * s;
    }
    template <template <class> class C, typename T>
    inline C<T> Abs(Tuple3<C, T> t) {
        return { std::abs(t.x), std::abs(t.y), std::abs(t.z) };
    }
    template <template <class> class C, typename T>
    inline C<T> Ceil(Tuple3<C, T> t) {
        return { std::ceil(t.x), std::ceil(t.y), std::ceil(t.z) };
    }
    template <template <class> class C, typename T>
    inline C<T> Floor(Tuple3<C, T> t) {
        return { std::floor(t.x), std::floor(t.y), std::floor(t.z) };
    }
    template <template <class> class C, typename T>
    inline auto Lerp(float t, Tuple3<C, T> t0, Tuple3<C, T> t1) {
        return (1 - t) * t0 + t * t1;
    }
    template <template <class> class C, typename T>
    inline C<T> FMA(float a, Tuple3<C, T> b, Tuple3<C, T> c) {
        return { std::fma(a, b.x, c.x), std::fma(a, b.y, c.y), std::fma(a, b.z, c.z) };
    }
    template <template <class> class C, typename T>
    inline C<T> FMA(Tuple3<C, T> a, float b, Tuple3<C, T> c) {
        return FMA(b, a, c);
    }
    template <template <class> class C, typename T>
    inline C<T> Min(Tuple3<C, T> t1, Tuple3<C, T> t2) {
        return { std::min(t1.x, t2.x), std::min(t1.y, t2.y), std::min(t1.z, t2.z) };
    }
    template <template <class> class C, typename T>
    inline T MinComponentValue(Tuple3<C, T> t) {
        return std::min({ t.x, t.y, t.z });
    }
    template <template <class> class C, typename T>
    inline int MinComponentIndex(Tuple3<C, T> t) {
        return (t.x < t.y) ? ((t.x < t.z) ? 0 : 2) : ((t.y < t.z) ? 1 : 2);
    }
    template <template <class> class C, typename T>
    inline C<T> Max(Tuple3<C, T> t1, Tuple3<C, T> t2) {
        return { std::max(t1.x, t2.x), std::max(t1.y, t2.y), std::max(t1.z, t2.z) };
    }
    template <template <class> class C, typename T>
    inline T MaxComponentValue(Tuple3<C, T> t) {
        return std::max({ t.x, t.y, t.z });
    }
    template <template <class> class C, typename T>
    inline int MaxComponentIndex(Tuple3<C, T> t) {
        return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2) : ((t.y > t.z) ? 1 : 2);
    }
    template <template <class> class C, typename T>
    inline C<T> Permute(Tuple3<C, T> t, std::array<int, 3> p) {
        return { t[p[0]], t[p[1]], t[p[2]] };
    }
    template <template <class> class C, typename T>
    inline T HProd(Tuple3<C, T> t) {
        return t.x * t.y * t.z;
    }

    // GeoVector2 Inline Functions
    template <typename T>
    inline auto Dot(GeoVector2<T> v1, GeoVector2<T> v2) -> typename TupleLength<T>::type {
        return SumOfProducts(v1.x, v2.x, v1.y, v2.y);
    }
    template <typename T>
    inline auto AbsDot(GeoVector2<T> v1, GeoVector2<T> v2) -> typename TupleLength<T>::type {
        return std::abs(Dot(v1, v2));
    }
    template <typename T>
    inline auto LengthSquared(GeoVector2<T> v) -> typename TupleLength<T>::type {
        return v.x * v.x + v.y * v.y;
    }
    template <typename T>
    inline auto Length(GeoVector2<T> v) -> typename TupleLength<T>::type {
        return std::sqrt(LengthSquared(v));
    }
    template <typename T>
    inline auto Normalize(GeoVector2<T> v) {
        return v / Length(v);
    }

    // GeoVector3 Inline Functions
    template <typename T>
    inline auto Dot(GeoVector3<T> v1, GeoVector3<T> v2) -> typename TupleLength<T>::type {
        return std::fma(v1.x, v2.x, SumOfProducts(v1.y, v2.y, v1.z, v2.z));
    }
    template <typename T>
    inline auto AbsDot(GeoVector3<T> v1, GeoVector3<T> v2) -> typename TupleLength<T>::type {
        return std::abs(Dot(v1, v2));
    }
    template <typename T>
    inline auto LengthSquared(GeoVector3<T> n) -> typename TupleLength<T>::type {
        return n.x * n.x + n.y * n.y + n.z * n.z;
    }
    template <typename T>
    inline auto Length(GeoVector3<T> v) -> typename TupleLength<T>::type {
        return std::sqrt(LengthSquared(v));
    }
    template <typename T>
    inline auto Normalize(GeoVector3<T> v) {
        return v / Length(v);
    }
    template <typename T>
    inline GeoVector3<T> Cross(GeoVector3<T> v1, GeoVector3<T> v2) {
        return { DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
                DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
                DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x) };
    }
    template <typename T>
    inline GeoVector3<T> FaceForward(GeoVector3<T> v1, GeoVector3<T> v2) {
        return (Dot(v1, v2) < 0.f) ? -v1 : v1;
    }
    // Equivalent to std::acos(Dot(a, b)), but more numerically stable.
    // via http://www.plunk.org/~hatch/rightway.html
    template <typename T>
    inline float AngleBetween(GeoVector3<T> v1, GeoVector3<T> v2) {
        if (Dot(v1, v2) < 0)
            return Pi - 2 * SafeASin(Length(v1 + v2) / 2);
        else
            return 2 * SafeASin(Length(v2 - v1) / 2);
    }

    template <typename T>
    inline TanVector3<T> GramSchmidt(TanVector3<T> v, TanVector3<T> w) {
        return v - Dot(v, w) * w;
    }
    template <typename T>
    inline void CoordinateSystem(GeoVector3<T> v1, TanVector3<T>* v2, TanVector3<T>* v3) {
        float sign = std::copysign(float(1), v1.z);
        float a = -1 / (sign + v1.z);
        float b = v1.x * v1.y * a;
        *v2 = TanVector3<T>(1 + sign * v1.x * v1.x * a, sign * b, -sign * v1.x);
        *v3 = TanVector3<T>(b, sign + v1.y * v1.y * a, -v1.y);
    }

    // Point Inline Functions
    template <typename T>
    inline auto Distance(Point2<T> p1, Point2<T> p2) -> typename TupleLength<T>::type {
        return Length(p1 - p2);
    }
    template <typename T>
    inline auto DistanceSquared(Point2<T> p1, Point2<T> p2) -> typename TupleLength<T>::type {
        return LengthSquared(p1 - p2);
    }
    template <typename T>
    inline auto Distance(Point3<T> p1, Point3<T> p2) {
        return Length(p1 - p2);
    }
    template <typename T>
    inline auto DistanceSquared(Point3<T> p1, Point3<T> p2) {
        return LengthSquared(p1 - p2);
    }

    // Quaternion Inline Functions
    inline Quaternion operator*(float f, Quaternion q) {
        return q * f;
    }
    inline float Dot(Quaternion q1, Quaternion q2) {
        return Dot(q1.v, q2.v) + q1.w * q2.w;
    }
    inline float Length(Quaternion q) {
        return std::sqrt(Dot(q, q));
    }
    inline Quaternion Normalize(Quaternion q) {
        return q / Length(q);
    }
    inline float AngleBetween(Quaternion q1, Quaternion q2) {
        if (Dot(q1, q2) < 0)
            return Pi - 2 * SafeASin(Length(q1 + q2) / 2);
        else
            return 2 * SafeASin(Length(q2 - q1) / 2);
    }
    // http://www.plunk.org/~hatch/rightway.html
    inline Quaternion Slerp(float t, Quaternion q1, Quaternion q2) {
        float theta = AngleBetween(q1, q2);
        float sinThetaOverTheta = SinXOverX(theta);
        return q1 * (1 - t) * SinXOverX((1 - t) * theta) / sinThetaOverTheta +
            q2 * t * SinXOverX(t * theta) / sinThetaOverTheta;
    }

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
        ret.pMin = b.pMin - GeoVector2<T>(delta, delta);
        ret.pMax = b.pMax + GeoVector2<T>(delta, delta);
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
        ret.pMin = b.pMin - Vector3<T>(delta, delta, delta);
        ret.pMax = b.pMax + Vector3<T>(delta, delta, delta);
        return ret;
    }
    template <typename T>
    inline bool Bounds3<T>::IntersectP(Point3f o, TanVector3f d, float tMax,
        float* hitt0, float* hitt1) const {
        float t0 = 0, t1 = tMax;
        for (int i = 0; i < 3; ++i) {
            // Update interval for _i_th bounding box slab
            float invRayDir = 1 / d[i];
            float tNear = (pMin[i] - o[i]) * invRayDir;
            float tFar = (pMax[i] - o[i]) * invRayDir;
            // Update parametric interval from slab intersection $t$ values
            if (tNear > tFar)
                pstd::swap(tNear, tFar);
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
    inline bool Bounds3<T>::IntersectP(Point3f o, TanVector3f d, float raytMax,
        TanVector3f invDir,
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
        // CotVectorly, the ending point is at the minimum x value and one past
        // the last valid y value.
        Point2i pEnd(b.pMin.x, b.pMax.y);
        // However, if the bounds are degenerate, override the end point to
        // equal the start point so that any attempt to iterate over the bounds
        // exits out immediately.
        if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
            pEnd = b.pMin;
        return Bounds2iIterator(b, pEnd);
    }

    // Spherical Geometry Inline Functions
    inline float SphericalTriangleArea(TanVector3f a, TanVector3f b, TanVector3f c) {
        return std::abs(
            2 * std::atan2(Dot(a, Cross(b, c)), 1 + Dot(a, b) + Dot(a, c) + Dot(b, c)));
    }
    inline float SphericalQuadArea(TanVector3f a, TanVector3f b, TanVector3f c, TanVector3f d) {
        TanVector3f axb = TanVector3f(Cross(a, b)), bxc = TanVector3f(Cross(b, c)),
            cxd = TanVector3f(Cross(c, d)), dxa = TanVector3f(Cross(d, a));
        if (LengthSquared(axb) == 0 || LengthSquared(bxc) == 0 || LengthSquared(cxd) == 0 ||
            LengthSquared(dxa) == 0)
            return 0;
        axb = TanVector3f(Normalize(axb));
        bxc = TanVector3f(Normalize(bxc));
        cxd = TanVector3f(Normalize(cxd));
        dxa = TanVector3f(Normalize(dxa));

        float alpha = AngleBetween(dxa, -axb);
        float beta = AngleBetween(axb, -bxc);
        float gamma = AngleBetween(bxc, -cxd);
        float delta = AngleBetween(cxd, -dxa);

        return std::abs(alpha + beta + gamma + delta - 2 * Pi);
    }
    inline GeoVector3f SphericalDirection(float sinTheta, float cosTheta, float phi) {
        return GeoVector3f(Clamp(sinTheta, -1, 1) * std::cos(phi),
            Clamp(sinTheta, -1, 1) * std::sin(phi), Clamp(cosTheta, -1, 1));
    }
    inline float SphericalTheta(GeoVector3f v) {
        return SafeACos(v.z);
    }
    inline float SphericalPhi(GeoVector3f v) {
        float p = std::atan2(v.y, v.x);
        return (p < 0) ? (p + 2 * Pi) : p;
    }
    inline float CosTheta(GeoVector3f w) {
        return w.z;
    }
    inline float Cos2Theta(GeoVector3f w) {
        return w.z * w.z;
    }
    inline float AbsCosTheta(GeoVector3f w) {
        return std::abs(w.z);
    }
    inline float Sin2Theta(GeoVector3f w) {
        return std::max<float>(0, 1 - Cos2Theta(w));
    }
    inline float SinTheta(GeoVector3f w) {
        return std::sqrt(Sin2Theta(w));
    }
    inline float TanTheta(GeoVector3f w) {
        return SinTheta(w) / CosTheta(w);
    }
    inline float Tan2Theta(GeoVector3f w) {
        return Sin2Theta(w) / Cos2Theta(w);
    }
    inline float CosPhi(GeoVector3f w) {
        float sinTheta = SinTheta(w);
        return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
    }
    inline float SinPhi(GeoVector3f w) {
        float sinTheta = SinTheta(w);
        return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
    }
    inline float CosDPhi(GeoVector3f wa, GeoVector3f wb) {
        float waxy = wa.x * wa.x + wa.y * wa.y, wbxy = wb.x * wb.x + wb.y * wb.y;
        if (waxy == 0 || wbxy == 0)
            return 1;
        return Clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt(waxy * wbxy), -1, 1);
    }
    inline bool SameHemisphere(TanVector3f w, TanVector3f wp) {
        return w.z * wp.z > 0;
    }
    inline bool SameHemisphere(TanVector3f w, CotVector3f wp) {
        return w.z * wp.z > 0;
    }

    class OctahedralVector {
    public:
        // OctahedralVector Public Methods
        OctahedralVector() = default;
        OctahedralVector(GeoVector3f v) {
            v /= std::abs(v.x) + std::abs(v.y) + std::abs(v.z);
            if (v.z >= 0) {
                x = Encode(v.x);
                y = Encode(v.y);
            }
            else {
                // Encode octahedral vector with $z < 0$
                x = Encode((1 - std::abs(v.y)) * Sign(v.x));
                y = Encode((1 - std::abs(v.x)) * Sign(v.y));
            }
        }

        explicit operator GeoVector3f() const {
            GeoVector3f v;
            v.x = -1 + 2 * (x / 65535.f);
            v.y = -1 + 2 * (y / 65535.f);
            v.z = 1 - (std::abs(v.x) + std::abs(v.y));
            // Reparameterize directions in the $z<0$ portion of the octahedron
            if (v.z < 0) {
                float xo = v.x;
                v.x = (1 - std::abs(v.y)) * Sign(xo);
                v.y = (1 - std::abs(xo)) * Sign(v.y);
            }

            return Normalize(v);
        }

    private:
        // OctahedralVector Private Methods
        static float Sign(float v) {
            return std::copysign(1.f, v);
        }

        static uint16_t Encode(float f) {
            return (uint16_t)std::round(Clamp((f + 1) / 2, 0, 1) * 65535.f);
        }

        // OctahedralVector Private Members
        uint16_t x, y;
    };

    class DirectionCone {
    public:
        // DirectionCone Public Methods
        DirectionCone() = default;
        DirectionCone(TanVector3f w, float cosTheta) : w(Normalize(w)), cosTheta(cosTheta) {}
        explicit DirectionCone(TanVector3f w) : DirectionCone(w, 1) {}

        bool IsEmpty() const {
            return cosTheta == Infinity;
        }
        static DirectionCone EntireSphere() {
            return DirectionCone(TanVector3f(0, 0, 1), -1);
        }

        TanVector3f ClosestVectorInCone(TanVector3f wp) const;

        // DirectionCone Public Members
        TanVector3f w;
        float cosTheta = Infinity;
    };

    // DirectionCone Function Declarations
    DirectionCone Union(const DirectionCone& a, const DirectionCone& b);

    // DirectionCone Inline Functions
    inline bool Inside(const DirectionCone& d, TanVector3f w) {
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
        TanVector3f w = TanVector3f(Normalize(pCenter - p));
        float sin2ThetaMax = radius * radius / DistanceSquared(pCenter, p);
        float cosThetaMax = SafeSqrt(1 - sin2ThetaMax);
        return DirectionCone(w, cosThetaMax);
    }
    inline TanVector3f DirectionCone::ClosestVectorInCone(TanVector3f wp) const {
        wp = TanVector3f(Normalize(wp));
        // Return provided vector if it is inside the cone
        if (Dot(wp, w) > cosTheta)
            return wp;

        // Find closest vector by rotating _wp_ until it touches the cone
        float sinTheta = -SafeSqrt(1 - cosTheta * cosTheta);
        TanVector3f a = TanVector3f(Cross(wp, w));
        return TanVector3f(cosTheta * w +
            (sinTheta / Length(a)) *
            TanVector3f(w.x * (wp.y * w.y + wp.z * w.z) - wp.x * (w.y * w.y + w.z * w.z),
                w.y * (wp.x * w.x + wp.z * w.z) - wp.y * (w.x * w.x + w.z * w.z),
                w.z * (wp.x * w.x + wp.y * w.y) - wp.z * (w.x * w.x + w.y * w.y)));
    }
    /*
    // https://www.iquilezles.org/www/articles/ibilinear/ibilinear.htm,
    // with a fix for perfect quads
    inline Point2f InvertBilinear(Point2f p, pstd::span<const Point2f> vert) {
        // The below assumes a quad (vs uv parametric layout) in v....
        Point2f a = vert[0], b = vert[1], c = vert[3], d = vert[2];
        GeoVector2f e = b - a, f = d - a, g = (a - b) + (c - d), h = p - a;

        auto cross2d = [](GeoVector2f a, GeoVector2f b) {
            return DifferenceOfProducts(a.x, b.y, a.y, b.x);
        };

        float k2 = cross2d(g, f);
        float k1 = cross2d(e, f) + cross2d(h, g);
        float k0 = cross2d(h, e);

        // if edges are parallel, this is a linear equation
        if (std::abs(k2) < 0.001f) {
            if (std::abs(e.x * k1 - g.x * k0) < 1e-5f)
                return Point2f((h.y * k1 + f.y * k0) / (e.y * k1 - g.y * k0), -k0 / k1);
            else
                return Point2f((h.x * k1 + f.x * k0) / (e.x * k1 - g.x * k0), -k0 / k1);
        }

        float v0, v1;
        if (!Quadratic(k2, k1, k0, &v0, &v1))
            return Point2f(0, 0);

        float u = (h.x - f.x * v0) / (e.x + g.x * v0);
        if (u < 0 || u > 1 || v0 < 0 || v0 > 1)
            return Point2f((h.x - f.x * v1) / (e.x + g.x * v1), v1);
        return Point2f(u, v0);
    }*/

    class Frame {
    public:
        // Frame Public Methods
        Frame() : x(1, 0, 0), y(0, 1, 0), z(0, 0, 1) {}
        Frame(TanVector3f x, TanVector3f y, TanVector3f z) : x(x), y(y), z(z) {}

        static Frame FromXZ(TanVector3f x, TanVector3f z) {
            return Frame(x, TanVector3f(Cross(z, x)), z);
        }
        static Frame FromXY(TanVector3f x, TanVector3f y) {
            return Frame(x, y, TanVector3f(Cross(x, y)));
        }
        static Frame FromZ(TanVector3f z) {
            TanVector3f x, y;
            CoordinateSystem(z, &x, &y);
            return Frame(x, y, z);
        }
        static Frame FromX(TanVector3f x) {
            TanVector3f y, z;
            CoordinateSystem(x, &y, &z);
            return Frame(x, y, z);
        }
        static Frame FromY(TanVector3f y) {
            TanVector3f x, z;
            CoordinateSystem(y, &z, &x);
            return Frame(x, y, z);
        }

        static Frame FromX(CotVector3f x) {
            TanVector3f y, z;
            CoordinateSystem(x, &y, &z);
            return Frame(TanVector3f(x), y, z);
        }
        static Frame FromY(CotVector3f y) {
            TanVector3f x, z;
            CoordinateSystem(y, &z, &x);
            return Frame(x, TanVector3f(y), z);
        }
        static Frame FromZ(CotVector3f z) {
            return FromZ(TanVector3f(z));
        }

        TanVector3f ToLocal(TanVector3f v) const {
            return TanVector3f(Dot(v, x), Dot(v, y), Dot(v, z));
        }
        CotVector3f ToLocal(CotVector3f n) const {
            return CotVector3f(Dot(n, x), Dot(n, y), Dot(n, z));
        }
        TanVector3f FromLocal(TanVector3f v) const {
            return TanVector3f(v.x * x + v.y * y + v.z * z);
        }
        CotVector3f FromLocal(CotVector3f n) const {
            return CotVector3f(n.x * x + n.y * y + n.z * z);
        }

        // Frame Public Members
        TanVector3f x, y, z;
    };

} //namespace lightfold