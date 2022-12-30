#pragma once

#include <cstdlib>
#include <cmath>
#include <array>

#include <math/math.h>

namespace lightfold {

namespace {
    // TupleLength Definition
    template <typename T>
    struct TupleLength {
        using type = float;
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

// Tuple2 Definition
template <template <typename> class Rep, typename T>
class Tuple2 {
public:
    static const int nDimensions = 2;

    // Tuple2 Public Methods
    // constructors
    Tuple2() = default;
    Tuple2(T x, T y) : x(x), y(y) {}
    
    // arithmetic operators
    template <typename U>
    auto operator+(Rep<U> c) const->Rep<decltype(T{} + U{})> {
        return { x + c.x, y + c.y };
    }
    template <typename U>
     Rep<T>& operator+=(Rep<U> c) {
        x += c.x;
        y += c.y;
        return static_cast<Rep<T> &>(*this);
    }
    template <typename U>
    auto operator-(Rep<U> c) const->Rep<decltype(T{} - U{})> {
        return { x - c.x, y - c.y };
    }
    template <typename U>
    Rep<T>& operator-=(Rep<U> c) {
        x -= c.x;
        y -= c.y;
        return static_cast<Rep<T> &>(*this);
    }
    template <typename U>
    auto operator*(U s) const->Rep<decltype(T{} *U{}) > {
        return { s * x, s * y };
    }
    template <typename U>
    Rep<T>& operator*=(U s) {
        x *= s;
        y *= s;
        return static_cast<Rep<T> &>(*this);
    }
    template <typename U>
    auto operator/(U d) const->Rep<decltype(T{} / U{}) > {
        return { x / d, y / d };
    }
    template <typename U>
    Rep<T>& operator/=(U d) {
        x /= d;
        y /= d;
        return static_cast<Rep<T> &>(*this);
    }
    Rep<T> operator-() const { return { -x, -y }; }

    // compare
    bool operator==(Rep<T> c) const { return x == c.x && y == c.y; }
    bool operator!=(Rep<T> c) const { return x != c.x || y != c.y; }

    // call members as in list
    T operator[](int i) const {
        return (i == 0) ? x : y;
    }
    T& operator[](int i) {
        return (i == 0) ? x : y;
    }

    // Tuple2 Public Members
    T x{}, y{};
};

// Tuple2 Inline Functions
template <template <class> class C, typename T, typename U>
inline auto operator*(U s, Tuple2<C, T> t)->C<decltype(T{} *U{}) > {
    return t * s;
}
template <template <class> class C, typename T>
inline C<T> Abs(Tuple2<C, T> t) {
    // "argument-dependent lookup..." (here and elsewhere)
    return { std::abs(t.x), std::abs(t.y) };
}
template <template <class> class C, typename T>
inline C<T> Ceil(Tuple2<C, T> t) {
    return { std::ceil(t.x), std::ceil(t.y) };
}
template <template <class> class C, typename T>
inline C<T> Floor(Tuple2<C, T> t) {
    return { std::floor(t.x), std::floor(t.y) };
}
template <template <class> class C, typename T>
inline auto Lerp(float t, Tuple2<C, T> t0, Tuple2<C, T> t1) {
    return (1 - t) * t0 + t * t1;
}
template <template <class> class C, typename T>
inline C<T> FMA(float a, Tuple2<C, T> b, Tuple2<C, T> c) {
    return { std::fma(a, b.x, c.x), std::fma(a, b.y, c.y) };
}
template <template <class> class C, typename T>
inline C<T> FMA(Tuple2<C, T> a, float b, Tuple2<C, T> c) {
    return FMA(b, a, c);
}
template <template <class> class C, typename T>
inline C<T> Min(Tuple2<C, T> t0, Tuple2<C, T> t1) {
    return { std::min(t0.x, t1.x), std::min(t0.y, t1.y) };
}
template <template <class> class C, typename T>
inline T MinComponentValue(Tuple2<C, T> t) {
    return std::min({ t.x, t.y });
}
template <template <class> class C, typename T>
inline int MinComponentIndex(Tuple2<C, T> t) {
    return (t.x < t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Max(Tuple2<C, T> t0, Tuple2<C, T> t1) {
    return { std::max(t0.x, t1.x), std::max(t0.y, t1.y) };
}

template <template <class> class C, typename T>
inline T MaxComponentValue(Tuple2<C, T> t) {
    return std::max({ t.x, t.y });
}

template <template <class> class C, typename T>
inline int MaxComponentIndex(Tuple2<C, T> t) {
    return (t.x > t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Permute(Tuple2<C, T> t, std::array<int, 2> p) {
    return { t[p[0]], t[p[1]] };
}

template <template <class> class C, typename T>
inline T HProd(Tuple2<C, T> t) {
    return t.x * t.y;
}

// Vector2 Definition
template <typename T>
class Vector2 : public Tuple2<Vector2, T> {
public:
    // Vector2 Public Methods
    using Tuple2<Vector2, T>::x;
    using Tuple2<Vector2, T>::y;

    Vector2() = default;
    Vector2(T x, T y) : Tuple2<pbrt::Vector2, T>(x, y) {}

    template <typename U>
    explicit Vector2(Vector2<U> v) : Tuple2<pbrt::Vector2, T>(T(v.x), T(v.y)) {}
};

// Vector2* Definitions
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;

// Vector2 Inline Functions
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

} //namespace lightfold