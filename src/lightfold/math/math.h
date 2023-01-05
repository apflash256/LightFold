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
        bool IsIdentity() const;

        std::span<const float> operator[](int i) const { return m[i]; }
        std::span<float> operator[](int i) { return std::span<float>(m[i]); }

    private:
        float m[N][N];
    };

    // SquareMatrix Inline Methods
    template <int N>
    inline bool SquareMatrix<N>::IsIdentity() const {
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

} // namespace lightfold