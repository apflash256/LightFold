#pragma once

#include <cstdlib>
#include <cmath>
#include <optional>
#include <span>

namespace lightfold {

    // Mathematical Constants
    constexpr float ShadowEpsilon = 0.0001f;
    constexpr float Pi = 3.14159265358979323846f;
    constexpr float InvPi = 0.31830988618379067154f;
    constexpr float Inv2Pi = 0.15915494309189533577f;
    constexpr float Inv4Pi = 0.07957747154594766788f;
    constexpr float PiOver2 = 1.57079632679489661923f;
    constexpr float PiOver4 = 0.78539816339744830961f;
    constexpr float Sqrt2 = 1.41421356237309504880f;

    static constexpr float Infinity = std::numeric_limits<float>::infinity();
    static constexpr float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5f;

    template <typename Ta, typename Tb, typename Tc, typename Td>
    inline auto DifferenceOfProducts(Ta a, Tb b, Tc c, Td d) {
        auto cd = c * d;
        auto differenceOfProducts = std::fma(a, b, -cd);
        auto error = std::fma(-c, d, cd);
        return differenceOfProducts + error;
    }

    template <typename Ta, typename Tb, typename Tc, typename Td>
    inline auto SumOfProducts(Ta a, Tb b, Tc c, Td d) {
        auto cd = c * d;
        auto sumOfProducts = std::fma(a, b, cd);
        auto error = std::fma(c, d, -cd);
        return sumOfProducts + error;
    }

    template <typename T, typename U, typename V>
    inline constexpr T Clamp(T val, U low, V high) {
        if (val < low)
            return T(low);
        else if (val > high)
            return T(high);
        else
            return val;
    }
    // http://www.plunk.org/~hatch/rightway.html
    inline float SinXOverX(float x) {
        if (1 + x * x == 1)
            return 1;
        return std::sin(x) / x;
    }
    inline float SafeASin(float x) {
        return std::asin(Clamp(x, -1, 1));
    }
    inline float SafeACos(float x) {
        return std::acos(Clamp(x, -1, 1));
    }
    inline double SafeASin(double x) {
        return std::asin(Clamp(x, -1, 1));
    }
    inline double SafeACos(double x) {
        return std::acos(Clamp(x, -1, 1));
    }

    inline float SafeSqrt(float x) {
        return std::sqrt(std::max(0.f, x));
    }
    inline double SafeSqrt(double x) {
        return std::sqrt(std::max(0., x));
    }

    inline constexpr float gamma(int n) {
        return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
    }

    namespace {
        template <int N>
        inline void init(float m[N][N], int i, int j) {}

        template <int N, typename... Args>
        inline void init(float m[N][N], int i, int j, float v, Args... args) {
            m[i][j] = v;
            if (++j == N) {
                ++i;
                j = 0;
            }
            init<N>(m, i, j, args...);
        }

        template <int N>
        inline void initDiag(float m[N][N], int i) {}

        template <int N, typename... Args>
        inline void initDiag(float m[N][N], int i, float v, Args... args) {
            m[i][i] = v;
            initDiag<N>(m, i + 1, args...);
        }
    } // anonymous namespace

    // SquareMatrix Definition
    template <int N>
    class SquareMatrix {
    public:
        // SquareMatrix Public Methods
        static SquareMatrix Zero() {
            SquareMatrix m;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    m.m[i][j] = 0;
            return m;
        }
        SquareMatrix() {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    m[i][j] = (i == j) ? 1.f : 0.f;
        }
        SquareMatrix(const float mat[N][N]) {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    m[i][j] = mat[i][j];
        }
        SquareMatrix(std::span<const float> t);

        template <typename... Args>
        SquareMatrix(float v, Args... args) {
            static_assert(1 + sizeof...(Args) == N * N,
                "Incorrect number of values provided to SquareMatrix constructor");
            init<N>(m, 0, 0, v, args...);
        }
        template <typename... Args>
        static SquareMatrix Diag(float v, Args... args) {
            static_assert(1 + sizeof...(Args) == N,
                "Incorrect number of values provided to SquareMatrix::Diag");
            SquareMatrix m;
            initDiag<N>(m.m, 0, v, args...);
            return m;
        }

        SquareMatrix operator+(const SquareMatrix& m) const {
            SquareMatrix r = *this;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    r.m[i][j] += m.m[i][j];
            return r;
        }
        SquareMatrix operator*(float s) const {
            SquareMatrix r = *this;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    r.m[i][j] *= s;
            return r;
        }
        SquareMatrix operator/(float s) const {
            SquareMatrix r = *this;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    r.m[i][j] /= s;
            return r;
        }
        bool operator==(const SquareMatrix<N>& m2) const {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    if (m[i][j] != m2.m[i][j])
                        return false;
            return true;
        }
        bool operator!=(const SquareMatrix<N>& m2) const {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    if (m[i][j] != m2.m[i][j])
                        return true;
            return false;
        }
        bool operator<(const SquareMatrix<N>& m2) const {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (m[i][j] < m2.m[i][j])
                        return true;
                    if (m[i][j] > m2.m[i][j])
                        return false;
                }
            return false;
        }
        bool IsIdentity() const {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) {
                        if (m[i][j] != 1)
                            return false;
                    }
                    else if (m[i][j] != 0)
                        return false;
                }
            return true;
        }

        std::span<const float> operator[](int i) const { return m[i]; }
        std::span<float> operator[](int i) { return std::span<float>(m[i]); }

    private:
        float m[N][N];
    };

    // SquareMatrix Inline Functions
    template <int N>
    inline SquareMatrix<N> operator*(float s, const SquareMatrix<N>& m) {
        return m * s;
    }

    template <typename Tresult, int N, typename T>
    inline Tresult Mul(const SquareMatrix<N>& m, const T& v) {
        Tresult result;
        for (int i = 0; i < N; ++i) {
            result[i] = 0;
            for (int j = 0; j < N; ++j)
                result[i] += m[i][j] * v[j];
        }
        return result;
    }

    template <int N>
    float Determinant(const SquareMatrix<N>& m);

    template <>
    inline float Determinant(const SquareMatrix<3>& m) {
        float minor12 = DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
        float minor02 = DifferenceOfProducts(m[1][0], m[2][2], m[1][2], m[2][0]);
        float minor01 = DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
        return std::fma(m[0][2], minor01,
            DifferenceOfProducts(m[0][0], minor12, m[0][1], minor02));
    }
    template <int N>
    SquareMatrix<N> InvertOrExit(const SquareMatrix<N>& m) {
        std::optional<SquareMatrix<N>> inv = Inverse(m);
        return *inv;
    }
    template <int N>
    inline SquareMatrix<N> Transpose(const SquareMatrix<N>& m) {
        SquareMatrix<N> r;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r[i][j] = m[j][i];
        return r;
    }
    inline std::optional<SquareMatrix<3>> Inverse(const SquareMatrix<3>& m) {
        float det = Determinant(m);
        if (det == 0)
            return {};
        float invDet = 1 / det;

        SquareMatrix<3> r;

        r[0][0] = invDet * DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
        r[1][0] = invDet * DifferenceOfProducts(m[1][2], m[2][0], m[1][0], m[2][2]);
        r[2][0] = invDet * DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
        r[0][1] = invDet * DifferenceOfProducts(m[0][2], m[2][1], m[0][1], m[2][2]);
        r[1][1] = invDet * DifferenceOfProducts(m[0][0], m[2][2], m[0][2], m[2][0]);
        r[2][1] = invDet * DifferenceOfProducts(m[0][1], m[2][0], m[0][0], m[2][1]);
        r[0][2] = invDet * DifferenceOfProducts(m[0][1], m[1][2], m[0][2], m[1][1]);
        r[1][2] = invDet * DifferenceOfProducts(m[0][2], m[1][0], m[0][0], m[1][2]);
        r[2][2] = invDet * DifferenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);

        return r;
    }

    template <int N, typename T>
    inline T operator*(const SquareMatrix<N>& m, const T& v) {
        return Mul<T>(m, v);
    }

    template <>
    inline SquareMatrix<4> operator*(const SquareMatrix<4>& m1,
        const SquareMatrix<4>& m2) {
        SquareMatrix<4> r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r[i][j] = m1[i][0] * m2[0][j] + m1[i][1] * m2[1][j] + m1[i][2] * m2[2][j] + m1[i][3] * m2[3][j];
        return r;
    }

    template <>
    inline SquareMatrix<3> operator*(const SquareMatrix<3>& m1,
        const SquareMatrix<3>& m2) {
        SquareMatrix<3> r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                r[i][j] = m1[i][0] * m2[0][j] + m1[i][1] * m2[1][j] + m1[i][2] * m2[2][j];
        return r;
    }

    template <int N>
    inline SquareMatrix<N> operator*(const SquareMatrix<N>& m1,
        const SquareMatrix<N>& m2) {
        SquareMatrix<N> r;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                r[i][j] = 0;
                for (int k = 0; k < N; ++k)
                    r[i][j] = std::fma(m1[i][k], m2[k][j], r[i][j]);
            }
        return r;
    }

    template <int N>
    inline SquareMatrix<N>::SquareMatrix(std::span<const float> t) {
        for (int i = 0; i < N * N; ++i)
            m[i / N][i % N] = t[i];
    }

    template <int N>
    SquareMatrix<N> operator*(const SquareMatrix<N>& m1,
        const SquareMatrix<N>& m2);

    template <>
    inline float Determinant(const SquareMatrix<1>& m) {
        return m[0][0];
    }

    template <>
    inline float Determinant(const SquareMatrix<2>& m) {
        return DifferenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);
    }

    template <>
    inline float Determinant(const SquareMatrix<4>& m) {
        float s0 = DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
        float s1 = DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
        float s2 = DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

        float s3 = DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
        float s4 = DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
        float s5 = DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

        float c0 = DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
        float c1 = DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
        float c2 = DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

        float c3 = DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
        float c4 = DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
        float c5 = DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

        return (DifferenceOfProducts(s0, c5, s1, c4) + DifferenceOfProducts(s2, c3, -s3, c2) +
            DifferenceOfProducts(s5, c0, s4, c1));
    }

    template <int N>
    inline float Determinant(const SquareMatrix<N>& m) {
        SquareMatrix<N - 1> sub;
        float det = 0;
        // Inefficient, but we don't currently use N>4 anyway..
        for (int i = 0; i < N; ++i) {
            // Sub-matrix without row 0 and column i
            for (int j = 0; j < N - 1; ++j)
                for (int k = 0; k < N - 1; ++k)
                    sub[j][k] = m[j + 1][k < i ? k : k + 1];

            float sign = (i & 1) ? -1 : 1;
            det += sign * m[0][i] * Determinant(sub);
        }
        return det;
    }

    inline std::optional<SquareMatrix<4>> Inverse(const SquareMatrix<4>& m) {
        // Via: https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
        // (c) Google, Apache license.

        // For 4x4 do not compute the adjugate as the transpose of the cofactor
        // matrix, because this results in extra work. Several calculations can be
        // shared across the sub-determinants.
        //
        // This approach is explained in David Eberly's Geometric Tools book,
        // excerpted here:
        //   http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
        float s0 = DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
        float s1 = DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
        float s2 = DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

        float s3 = DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
        float s4 = DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
        float s5 = DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

        float c0 = DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
        float c1 = DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
        float c2 = DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

        float c3 = DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
        float c4 = DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
        float c5 = DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

        float determinant = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 + s5 * c0 - s4 * c1;
        if (determinant == 0)
            return {};
        float s = 1 / determinant;

        float inv[4][4] = { {s * (m[1][1] * c5 + m[1][3] * c3 - m[1][2] * c4),
                            s * (-m[0][1] * c5 + m[0][2] * c4 - m[0][3] * c3),
                            s * (m[3][1] * s5 + m[3][3] * s3 - m[3][2] * s4),
                            s * (-m[2][1] * s5 + m[2][2] * s4 - m[2][3] * s3)},

                           {s * (-m[1][0] * c5 + m[1][2] * c2 - m[1][3] * c1),
                            s * (m[0][0] * c5 + m[0][3] * c1 - m[0][2] * c2),
                            s * (-m[3][0] * s5 + m[3][2] * s2 - m[3][3] * s1),
                            s * (m[2][0] * s5 + m[2][3] * s1 - m[2][2] * s2)},

                           {s * (m[1][0] * c4 + m[1][3] * c0 - m[1][1] * c2),
                            s * (-m[0][0] * c4 + m[0][1] * c2 - m[0][3] * c0),
                            s * (m[3][0] * s4 + m[3][3] * s0 - m[3][1] * s2),
                            s * (-m[2][0] * s4 + m[2][1] * s2 - m[2][3] * s0)},

                           {s * (-m[1][0] * c3 + m[1][1] * c1 - m[1][2] * c0),
                            s * (m[0][0] * c3 + m[0][2] * c0 - m[0][1] * c1),
                            s * (-m[3][0] * s3 + m[3][1] * s1 - m[3][2] * s0),
                            s * (m[2][0] * s3 + m[2][2] * s0 - m[2][1] * s1)} };

        return SquareMatrix<4>(inv);
    }

    extern template class SquareMatrix<2>;
    extern template class SquareMatrix<3>;
    extern template class SquareMatrix<4>;

    // Interval Definition
    class Interval {
    public:
        // Interval Public Methods
        Interval() = default;
        explicit Interval(float v) : low(v), high(v) {}
        constexpr Interval(float low, float high) : low(std::min(low, high)), high(std::max(low, high)) {}

        static Interval FromValueAndError(float v, float err) {
            Interval i;
            if (err == 0)
                i.low = i.high = v;
            else {
                i.low = v - err;
                i.high = v + err;
            }
            return i;
        }
        Interval& operator=(float v) {
            low = high = v;
            return *this;
        }

        float UpperBound() const { return high; }
        float LowerBound() const { return low; }
        float Midpoint() const { return (low + high) / 2; }
        float Width() const { return high - low; }

        float operator[](int i) const {
            return (i == 0) ? low : high;
        }
        explicit operator float() const { return Midpoint(); }

        bool Exactly(float v) const { return low == v && high == v; }
        bool operator==(float v) const { return Exactly(v); }

        Interval operator-() const { return { -high, -low }; }
        Interval operator+(Interval i) const {
            return { low + i.low, high + i.high };
        }
        Interval operator-(Interval i) const {
            return { low - i.high, high - i.low };
        }
        Interval operator*(Interval i) const {
            float lp[4] = { low * i.low, high * i.low,
                           low * i.high, high * i.high };
            float hp[4] = { low * i.low, high * i.low,
                           low * i.high, high * i.high };
            return { std::min({lp[0], lp[1], lp[2], lp[3]}),
                    std::max({hp[0], hp[1], hp[2], hp[3]}) };
        }
        Interval operator/(Interval i) const;

        bool operator==(Interval i) const {
            return low == i.low && high == i.high;
        }
        bool operator!=(float f) const { return f < low || f > high; }

        Interval& operator+=(Interval i) {
            *this = Interval(*this + i);
            return *this;
        }
        Interval& operator-=(Interval i) {
            *this = Interval(*this - i);
            return *this;
        }
        Interval& operator*=(Interval i) {
            *this = Interval(*this * i);
            return *this;
        }
        Interval& operator/=(Interval i) {
            *this = Interval(*this / i);
            return *this;
        }

        Interval& operator+=(float f) { return *this += Interval(f); }
        Interval& operator-=(float f) { return *this -= Interval(f); }
        Interval& operator*=(float f) {
            if (f > 0)
                *this = Interval(f * low, f * high);
            else
                *this = Interval(f * high, f * low);
            return *this;
        }

        Interval& operator/=(float f) {
            if (f > 0)
                *this = Interval(low / f, high / f);
            else
                *this = Interval(high / f, low / f);
            return *this;
        }

        // Interval Private Members
    private:
        float low, high;
    };

    // Interval Inline Functions
    inline bool InRange(float v, Interval i) {
        return v >= i.LowerBound() && v <= i.UpperBound();
    }
    inline bool InRange(Interval a, Interval b) {
        return a.LowerBound() <= b.UpperBound() && a.UpperBound() >= b.LowerBound();
    }

    inline Interval Interval::operator/(Interval i) const {
        if (InRange(0, i))
            // The interval we're dividing by straddles zero, so just
            // return an interval of everything.
            return Interval(-Infinity, Infinity);

        float lowQuot[4] = { low / i.low, high / i.low, low / i.high, high / i.high };
        float highQuot[4] = { low / i.low, high / i.low, low / i.high, high / i.high };
        return { std::min({lowQuot[0], lowQuot[1], lowQuot[2], lowQuot[3]}),
                std::max({highQuot[0], highQuot[1], highQuot[2], highQuot[3]}) };
    }

    inline Interval Sqr(Interval i) {
        float alow = std::abs(i.LowerBound()), ahigh = std::abs(i.UpperBound());
        if (alow > ahigh)
            std::swap(alow, ahigh);
        if (InRange(0, i))
            return Interval(0, ahigh * ahigh);
        return Interval(alow * alow, ahigh * ahigh);
    }

    inline Interval MulPow2(float s, Interval i);
    inline Interval MulPow2(Interval i, float s);

    inline Interval operator+(float f, Interval i) {
        return Interval(f) + i;
    }

    inline Interval operator-(float f, Interval i) {
        return Interval(f) - i;
    }

    inline Interval operator*(float f, Interval i) {
        if (f > 0)
            return Interval(f * i.LowerBound(), f * i.UpperBound());
        else
            return Interval(f * i.UpperBound(), f * i.LowerBound());
    }

    inline Interval operator/(float f, Interval i) {
        if (InRange(0, i))
            // The interval we're dividing by straddles zero, so just
            // return an interval of everything.
            return Interval(-Infinity, Infinity);

        if (f > 0)
            return Interval(f / i.UpperBound(), f / i.LowerBound());
        else
            return Interval(f / i.LowerBound(), f / i.UpperBound());
    }

    inline Interval operator+(Interval i, float f) {
        return i + Interval(f);
    }

    inline Interval operator-(Interval i, float f) {
        return i - Interval(f);
    }

    inline Interval operator*(Interval i, float f) {
        if (f > 0)
            return Interval(f * i.LowerBound(), f * i.UpperBound());
        else
            return Interval(f * i.UpperBound(), f * i.LowerBound());
    }

    inline Interval operator/(Interval i, float f) {
        if (f == 0)
            return Interval(-Infinity, Infinity);

        if (f > 0)
            return Interval(i.LowerBound() / f, i.UpperBound() / f);
        else
            return Interval(i.UpperBound() / f, i.LowerBound() / f);
    }

    inline float Floor(Interval i) {
        return std::floor(i.LowerBound());
    }

    inline float Ceil(Interval i) {
        return std::ceil(i.UpperBound());
    }

    inline float floor(Interval i) {
        return Floor(i);
    }

    inline float ceil(Interval i) {
        return Ceil(i);
    }

    inline float Min(Interval a, Interval b) {
        return std::min(a.LowerBound(), b.LowerBound());
    }

    inline float Max(Interval a, Interval b) {
        return std::max(a.UpperBound(), b.UpperBound());
    }

    inline Interval Sqrt(Interval i) {
        return { std::sqrt(i.LowerBound()), std::sqrt(i.UpperBound()) };
    }

    inline Interval FMA(Interval a, Interval b, Interval c) {
        float low = std::min({ std::fma(a.LowerBound(), b.LowerBound(), c.LowerBound()),
                              std::fma(a.UpperBound(), b.LowerBound(), c.LowerBound()),
                              std::fma(a.LowerBound(), b.UpperBound(), c.LowerBound()),
                              std::fma(a.UpperBound(), b.UpperBound(), c.LowerBound()) });
        float high = std::max({ std::fma(a.LowerBound(), b.LowerBound(), c.UpperBound()),
                               std::fma(a.UpperBound(), b.LowerBound(), c.UpperBound()),
                               std::fma(a.LowerBound(), b.UpperBound(), c.UpperBound()),
                               std::fma(a.UpperBound(), b.UpperBound(), c.UpperBound()) });
        return Interval(low, high);
    }

    inline Interval DifferenceOfProducts(Interval a, Interval b, Interval c,
        Interval d) {
        float ab[4] = { a.LowerBound() * b.LowerBound(), a.UpperBound() * b.LowerBound(),
                       a.LowerBound() * b.UpperBound(), a.UpperBound() * b.UpperBound() };
        float abLow = std::min({ ab[0], ab[1], ab[2], ab[3] });
        float abHigh = std::max({ ab[0], ab[1], ab[2], ab[3] });
        int abLowIndex = abLow == ab[0] ? 0 : (abLow == ab[1] ? 1 : (abLow == ab[2] ? 2 : 3));
        int abHighIndex =
            abHigh == ab[0] ? 0 : (abHigh == ab[1] ? 1 : (abHigh == ab[2] ? 2 : 3));

        float cd[4] = { c.LowerBound() * d.LowerBound(), c.UpperBound() * d.LowerBound(),
                       c.LowerBound() * d.UpperBound(), c.UpperBound() * d.UpperBound() };
        float cdLow = std::min({ cd[0], cd[1], cd[2], cd[3] });
        float cdHigh = std::max({ cd[0], cd[1], cd[2], cd[3] });
        int cdLowIndex = cdLow == cd[0] ? 0 : (cdLow == cd[1] ? 1 : (cdLow == cd[2] ? 2 : 3));
        int cdHighIndex =
            cdHigh == cd[0] ? 0 : (cdHigh == cd[1] ? 1 : (cdHigh == cd[2] ? 2 : 3));

        // Invert cd Indices since it's subtracted...
        float low = DifferenceOfProducts(a[abLowIndex & 1], b[abLowIndex >> 1],
            c[cdHighIndex & 1], d[cdHighIndex >> 1]);
        float high = DifferenceOfProducts(a[abHighIndex & 1], b[abHighIndex >> 1],
            c[cdLowIndex & 1], d[cdLowIndex >> 1]);

        return { low, high };
    }

    inline Interval SumOfProducts(Interval a, Interval b, Interval c,
        Interval d) {
        return DifferenceOfProducts(a, b, -c, d);
    }

    inline Interval MulPow2(float s, Interval i) {
        return MulPow2(i, s);
    }

    inline Interval MulPow2(Interval i, float s) {
        float as = std::abs(s);
        // Multiplication by powers of 2 is exaact
        return Interval(std::min(i.LowerBound() * s, i.UpperBound() * s),
            std::max(i.LowerBound() * s, i.UpperBound() * s));
    }

    inline Interval Abs(Interval i) {
        if (i.LowerBound() >= 0)
            // The entire interval is greater than zero, so we're all set.
            return i;
        else if (i.UpperBound() <= 0)
            // The entire interval is less than zero.
            return Interval(-i.UpperBound(), -i.LowerBound());
        else
            // The interval straddles zero.
            return Interval(0, std::max(-i.LowerBound(), i.UpperBound()));
    }

    inline Interval abs(Interval i) {
        return Abs(i);
    }

    inline Interval ACos(Interval i) {
        float low = std::acos(std::min<float>(1, i.UpperBound()));
        float high = std::acos(std::max<float>(-1, i.LowerBound()));

        return Interval(std::max<float>(0, low), high);
    }

    inline Interval Sin(Interval i) {
        float low = std::sin(std::max<float>(0, i.LowerBound()));
        float high = std::sin(i.UpperBound());
        if (low > high)
            std::swap(low, high);
        low = std::max<float>(-1, low);
        high = std::min<float>(1, high);
        if (InRange(Pi / 2, i))
            high = 1;
        if (InRange((3.f / 2.f) * Pi, i))
            low = -1;

        return Interval(low, high);
    }

    inline Interval Cos(Interval i) {
        float low = std::cos(std::max<float>(0, i.LowerBound()));
        float high = std::cos(i.UpperBound());
        if (low > high)
            std::swap(low, high);
        low = std::max<float>(-1, low);
        high = std::min<float>(1, high);
        if (InRange(Pi, i))
            low = -1;

        return Interval(low, high);
    }

    inline bool Quadratic(Interval a, Interval b, Interval c, Interval* t0,
        Interval* t1) {
        // Find quadratic discriminant
        Interval discrim = DifferenceOfProducts(b, b, MulPow2(4, a), c);
        if (discrim.LowerBound() < 0)
            return false;
        Interval floatRootDiscrim = Sqrt(discrim);

        // Compute quadratic _t_ values
        Interval q;
        if ((float)b < 0)
            q = MulPow2(-.5, b - floatRootDiscrim);
        else
            q = MulPow2(-.5, b + floatRootDiscrim);
        *t0 = q / a;
        *t1 = c / q;
        if (t0->LowerBound() > t1->LowerBound())
            std::swap(*t0, *t1);
        return true;
    }

    inline Interval SumSquares(Interval i) {
        return Sqr(i);
    }

    template <typename... Args>
    inline Interval SumSquares(Interval i, Args... args) {
        Interval ss = FMA(i, i, SumSquares(args...));
        return Interval(std::max<float>(0, ss.LowerBound()), ss.UpperBound());
    }

    static constexpr int PrimeTableSize = 1000;
    static const int Primes[PrimeTableSize] = {
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
    193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
    397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
    499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
    727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829,
    839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953,
    967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051,
    1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259,
    1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367,
    1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471,
    1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567,
    1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
    1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777,
    1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879,
    1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099,
    2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221,
    2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333,
    2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417,
    2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549,
    2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671,
    2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
    2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857,
    2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971,
    2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109,
    3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229,
    3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343,
    3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463,
    3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559,
    3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673,
    3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911,
    3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019,
    4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133,
    4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253,
    4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373,
    4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507,
    4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637,
    4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877,
    4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987,
    4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099,
    5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231,
    5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381,
    5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477,
    5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581,
    5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827,
    5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927,
    5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079,
    6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203,
    6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311,
    6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421,
    6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569,
    6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691,
    6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823,
    6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947,
    6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039,
    7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193,
    7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321,
    7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481,
    7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573,
    7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
    7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823,
    7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 };

    static const int PrimeSums[PrimeTableSize] = {
    0, 2, 5, 10, 17, 28, 41, 58, 77, 100, 129, 160, 197, 238, 281, 328, 381, 440,
    501, 568, 639, 712, 791, 874, 963, 1060, 1161, 1264, 1371, 1480, 1593, 1720,
    1851, 1988, 2127, 2276, 2427, 2584, 2747, 2914, 3087, 3266, 3447, 3638,
    3831, 4028, 4227, 4438, 4661, 4888, 5117, 5350, 5589, 5830, 6081, 6338,
    6601, 6870, 7141, 7418, 7699, 7982, 8275, 8582, 8893, 9206, 9523, 9854,
    10191, 10538, 10887, 11240, 11599, 11966, 12339, 12718, 13101, 13490, 13887,
    14288, 14697, 15116, 15537, 15968, 16401, 16840, 17283, 17732, 18189, 18650,
    19113, 19580, 20059, 20546, 21037, 21536, 22039, 22548, 23069, 23592, 24133,
    24680, 25237, 25800, 26369, 26940, 27517, 28104, 28697, 29296, 29897, 30504,
    31117, 31734, 32353, 32984, 33625, 34268, 34915, 35568, 36227, 36888, 37561,
    38238, 38921, 39612, 40313, 41022, 41741, 42468, 43201, 43940, 44683, 45434,
    46191, 46952, 47721, 48494, 49281, 50078, 50887, 51698, 52519, 53342, 54169,
    54998, 55837, 56690, 57547, 58406, 59269, 60146, 61027, 61910, 62797, 63704,
    64615, 65534, 66463, 67400, 68341, 69288, 70241, 71208, 72179, 73156, 74139,
    75130, 76127, 77136, 78149, 79168, 80189, 81220, 82253, 83292, 84341, 85392,
    86453, 87516, 88585, 89672, 90763, 91856, 92953, 94056, 95165, 96282, 97405,
    98534, 99685, 100838, 102001, 103172, 104353, 105540, 106733, 107934,
    109147, 110364, 111587, 112816, 114047, 115284, 116533, 117792, 119069,
    120348, 121631, 122920, 124211, 125508, 126809, 128112, 129419, 130738,
    132059, 133386, 134747, 136114, 137487, 138868, 140267, 141676, 143099,
    144526, 145955, 147388, 148827, 150274, 151725, 153178, 154637, 156108,
    157589, 159072, 160559, 162048, 163541, 165040, 166551, 168074, 169605,
    171148, 172697, 174250, 175809, 177376, 178947, 180526, 182109, 183706,
    185307, 186914, 188523, 190136, 191755, 193376, 195003, 196640, 198297,
    199960, 201627, 203296, 204989, 206686, 208385, 210094, 211815, 213538,
    215271, 217012, 218759, 220512, 222271, 224048, 225831, 227618, 229407,
    231208, 233019, 234842, 236673, 238520, 240381, 242248, 244119, 245992,
    247869, 249748, 251637, 253538, 255445, 257358, 259289, 261222, 263171,
    265122, 267095, 269074, 271061, 273054, 275051, 277050, 279053, 281064,
    283081, 285108, 287137, 289176, 291229, 293292, 295361, 297442, 299525,
    301612, 303701, 305800, 307911, 310024, 312153, 314284, 316421, 318562,
    320705, 322858, 325019, 327198, 329401, 331608, 333821, 336042, 338279,
    340518, 342761, 345012, 347279, 349548, 351821, 354102, 356389, 358682,
    360979, 363288, 365599, 367932, 370271, 372612, 374959, 377310, 379667,
    382038, 384415, 386796, 389179, 391568, 393961, 396360, 398771, 401188,
    403611, 406048, 408489, 410936, 413395, 415862, 418335, 420812, 423315,
    425836, 428367, 430906, 433449, 435998, 438549, 441106, 443685, 446276,
    448869, 451478, 454095, 456716, 459349, 461996, 464653, 467312, 469975,
    472646, 475323, 478006, 480693, 483382, 486075, 488774, 491481, 494192,
    496905, 499624, 502353, 505084, 507825, 510574, 513327, 516094, 518871,
    521660, 524451, 527248, 530049, 532852, 535671, 538504, 541341, 544184,
    547035, 549892, 552753, 555632, 558519, 561416, 564319, 567228, 570145,
    573072, 576011, 578964, 581921, 584884, 587853, 590824, 593823, 596824,
    599835, 602854, 605877, 608914, 611955, 615004, 618065, 621132, 624211,
    627294, 630383, 633492, 636611, 639732, 642869, 646032, 649199, 652368,
    655549, 658736, 661927, 665130, 668339, 671556, 674777, 678006, 681257,
    684510, 687767, 691026, 694297, 697596, 700897, 704204, 707517, 710836,
    714159, 717488, 720819, 724162, 727509, 730868, 734229, 737600, 740973,
    744362, 747753, 751160, 754573, 758006, 761455, 764912, 768373, 771836,
    775303, 778772, 782263, 785762, 789273, 792790, 796317, 799846, 803379,
    806918, 810459, 814006, 817563, 821122, 824693, 828274, 831857, 835450,
    839057, 842670, 846287, 849910, 853541, 857178, 860821, 864480, 868151,
    871824, 875501, 879192, 882889, 886590, 890299, 894018, 897745, 901478,
    905217, 908978, 912745, 916514, 920293, 924086, 927883, 931686, 935507,
    939330, 943163, 947010, 950861, 954714, 958577, 962454, 966335, 970224,
    974131, 978042, 981959, 985878, 989801, 993730, 997661, 1001604, 1005551,
    1009518, 1013507, 1017508, 1021511, 1025518, 1029531, 1033550, 1037571,
    1041598, 1045647, 1049698, 1053755, 1057828, 1061907, 1065998, 1070091,
    1074190, 1078301, 1082428, 1086557, 1090690, 1094829, 1098982, 1103139,
    1107298, 1111475, 1115676, 1119887, 1124104, 1128323, 1132552, 1136783,
    1141024, 1145267, 1149520, 1153779, 1158040, 1162311, 1166584, 1170867,
    1175156, 1179453, 1183780, 1188117, 1192456, 1196805, 1201162, 1205525,
    1209898, 1214289, 1218686, 1223095, 1227516, 1231939, 1236380, 1240827,
    1245278, 1249735, 1254198, 1258679, 1263162, 1267655, 1272162, 1276675,
    1281192, 1285711, 1290234, 1294781, 1299330, 1303891, 1308458, 1313041,
    1317632, 1322229, 1326832, 1331453, 1336090, 1340729, 1345372, 1350021,
    1354672, 1359329, 1363992, 1368665, 1373344, 1378035, 1382738, 1387459,
    1392182, 1396911, 1401644, 1406395, 1411154, 1415937, 1420724, 1425513,
    1430306, 1435105, 1439906, 1444719, 1449536, 1454367, 1459228, 1464099,
    1468976, 1473865, 1478768, 1483677, 1488596, 1493527, 1498460, 1503397,
    1508340, 1513291, 1518248, 1523215, 1528184, 1533157, 1538144, 1543137,
    1548136, 1553139, 1558148, 1563159, 1568180, 1573203, 1578242, 1583293,
    1588352, 1593429, 1598510, 1603597, 1608696, 1613797, 1618904, 1624017,
    1629136, 1634283, 1639436, 1644603, 1649774, 1654953, 1660142, 1665339,
    1670548, 1675775, 1681006, 1686239, 1691476, 1696737, 1702010, 1707289,
    1712570, 1717867, 1723170, 1728479, 1733802, 1739135, 1744482, 1749833,
    1755214, 1760601, 1765994, 1771393, 1776800, 1782213, 1787630, 1793049,
    1798480, 1803917, 1809358, 1814801, 1820250, 1825721, 1831198, 1836677,
    1842160, 1847661, 1853164, 1858671, 1864190, 1869711, 1875238, 1880769,
    1886326, 1891889, 1897458, 1903031, 1908612, 1914203, 1919826, 1925465,
    1931106, 1936753, 1942404, 1948057, 1953714, 1959373, 1965042, 1970725,
    1976414, 1982107, 1987808, 1993519, 1999236, 2004973, 2010714, 2016457,
    2022206, 2027985, 2033768, 2039559, 2045360, 2051167, 2056980, 2062801,
    2068628, 2074467, 2080310, 2086159, 2092010, 2097867, 2103728, 2109595,
    2115464, 2121343, 2127224, 2133121, 2139024, 2144947, 2150874, 2156813,
    2162766, 2168747, 2174734, 2180741, 2186752, 2192781, 2198818, 2204861,
    2210908, 2216961, 2223028, 2229101, 2235180, 2241269, 2247360, 2253461,
    2259574, 2265695, 2271826, 2277959, 2284102, 2290253, 2296416, 2302589,
    2308786, 2314985, 2321188, 2327399, 2333616, 2339837, 2346066, 2352313,
    2358570, 2364833, 2371102, 2377373, 2383650, 2389937, 2396236, 2402537,
    2408848, 2415165, 2421488, 2427817, 2434154, 2440497, 2446850, 2453209,
    2459570, 2465937, 2472310, 2478689, 2485078, 2491475, 2497896, 2504323,
    2510772, 2517223, 2523692, 2530165, 2536646, 2543137, 2549658, 2556187,
    2562734, 2569285, 2575838, 2582401, 2588970, 2595541, 2602118, 2608699,
    2615298, 2621905, 2628524, 2635161, 2641814, 2648473, 2655134, 2661807,
    2668486, 2675175, 2681866, 2688567, 2695270, 2701979, 2708698, 2715431,
    2722168, 2728929, 2735692, 2742471, 2749252, 2756043, 2762836, 2769639,
    2776462, 2783289, 2790118, 2796951, 2803792, 2810649, 2817512, 2824381,
    2831252, 2838135, 2845034, 2851941, 2858852, 2865769, 2872716, 2879665,
    2886624, 2893585, 2900552, 2907523, 2914500, 2921483, 2928474, 2935471,
    2942472, 2949485, 2956504, 2963531, 2970570, 2977613, 2984670, 2991739,
    2998818, 3005921, 3013030, 3020151, 3027278, 3034407, 3041558, 3048717,
    3055894, 3063081, 3070274, 3077481, 3084692, 3091905, 3099124, 3106353,
    3113590, 3120833, 3128080, 3135333, 3142616, 3149913, 3157220, 3164529,
    3171850, 3179181, 3186514, 3193863, 3201214, 3208583, 3215976, 3223387,
    3230804, 3238237, 3245688, 3253145, 3260604, 3268081, 3275562, 3283049,
    3290538, 3298037, 3305544, 3313061, 3320584, 3328113, 3335650, 3343191,
    3350738, 3358287, 3365846, 3373407, 3380980, 3388557, 3396140, 3403729,
    3411320, 3418923, 3426530, 3434151, 3441790, 3449433, 3457082, 3464751,
    3472424, 3480105, 3487792, 3495483, 3503182, 3510885, 3518602, 3526325,
    3534052, 3541793, 3549546, 3557303, 3565062, 3572851, 3580644, 3588461,
    3596284, 3604113, 3611954, 3619807, 3627674, 3635547, 3643424, 3651303,
    3659186, 3667087, 3674994 };
} // namespace lightfold