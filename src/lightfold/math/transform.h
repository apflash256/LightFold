#pragma once

#include <math/ray.h>

namespace lightfold {

    class SurfaceInteraction;

    // Quaternion Definition
    class Quaternion {
    public:
        // Quaternion Public Methods
        Quaternion() = default;
        Quaternion(Tangent3f v, float w) : v(v), w(w) {}

        Quaternion& operator+=(Quaternion q) {
            v += q.v;
            w += q.w;
            return *this;
        }
        Quaternion& operator-=(Quaternion q) {
            v -= q.v;
            w -= q.w;
            return *this;
        }
        Quaternion& operator*=(float f) {
            v *= f;
            w *= f;
            return *this;
        }
        Quaternion& operator/=(float f) {
            v /= f;
            w /= f;
            return *this;
        }
        Quaternion operator+(Quaternion q) const {
            return { v + q.v, w + q.w };
        }
        Quaternion operator-(Quaternion q) const {
            return { v - q.v, w - q.w };
        }
        Quaternion operator*(float f) const {
            return { v * f, w * f };
        }
        Quaternion operator/(float f) const {
            return { v / f, w / f };
        }
        Quaternion operator-() const {
            return { -v, -w };
        }

        // Quaternion Public Members
        Tangent3f v;
        float w = 1;
    };

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

    // Transform Definition
    class Transform {
    public:
        // Transform Public Methods
        template <typename T>
        inline Tangent3<T> ApplyInverse(Tangent3<T> v) const;
        template <typename T>
        inline Normal3<T> ApplyInverse(Normal3<T> v) const;
        template <typename T>
        inline Point3<T> ApplyInverse(Point3<T> p) const;
        inline Ray ApplyInverse(const Ray& r, float* tMax = nullptr) const;
        inline RayDifferential ApplyInverse(const RayDifferential& r, float* tMax = nullptr) const;

        Transform() = default;
        Transform(const SquareMatrix<4>& m) : m(m) {
            std::optional<SquareMatrix<4>> inv = Inverse(m);
            if (inv)
                mInv = *inv;
            else {
                // Initialize _mInv_ with not-a-number values
                float NaN = std::numeric_limits<float>::has_signaling_NaN
                    ? std::numeric_limits<float>::signaling_NaN()
                    : std::numeric_limits<float>::quiet_NaN();
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        mInv[i][j] = NaN;
            }
        }
        Transform(const float mat[4][4]) : Transform(SquareMatrix<4>(mat)) {}
        Transform(const SquareMatrix<4>& m, const SquareMatrix<4>& mInv) : m(m), mInv(mInv) {}

        const SquareMatrix<4>& GetMatrix() const { return m; }
        const SquareMatrix<4>& GetInverseMatrix() const { return mInv; }

        bool operator==(const Transform& t) const { return t.m == m; }
        bool operator!=(const Transform& t) const { return t.m != m; }
        bool IsIdentity() const { return m.IsIdentity(); }
        bool HasScale(float tolerance = 1e-3f) const {
            float la2 = LengthSquared((*this)(Tangent3f(1, 0, 0)));
            float lb2 = LengthSquared((*this)(Tangent3f(0, 1, 0)));
            float lc2 = LengthSquared((*this)(Tangent3f(0, 0, 1)));
            return (std::abs(la2 - 1) > tolerance || std::abs(lb2 - 1) > tolerance ||
                std::abs(lc2 - 1) > tolerance);
        }

        template <typename T>
        Point3<T> operator()(Point3<T> p) const;
        template <typename T>
        Tangent3<T> operator()(Tangent3<T> v) const;
        template <typename T>
        Normal3<T> operator()(Normal3<T>) const;
        Bounds3f operator()(const Bounds3f& b) const;
        Ray operator()(const Ray& r, float* tMax = nullptr) const;
        RayDifferential operator()(const RayDifferential& r, float* tMax = nullptr) const;
        SurfaceInteraction operator()(const SurfaceInteraction& si) const;
        
        template <typename T>
        inline Point3<T> operator()(const Point3<T>& pt,
            Tangent3<T>* absError) const;
        template <typename T>
        inline Point3<T> operator()(const Point3<T>& p, const Tangent3<T>& pError,
            Tangent3<T>* pTransError) const;
        template <typename T>
        inline Tangent3<T> operator()(const Tangent3<T>& v,
            Tangent3<T>* vTransError) const;
        template <typename T>
        inline Tangent3<T> operator()(const Tangent3<T>& v, const Tangent3<T>& vError,
            Tangent3<T>* vTransError) const;
        inline Ray operator()(const Ray& r, Tangent3f* oError,
            Tangent3f* dError) const;
        inline Ray operator()(const Ray& r, const Tangent3f& oErrorIn,
            const Tangent3f& dErrorIn, Tangent3f* oErrorOut,
            Tangent3f* dErrorOut) const;

        Transform operator*(const Transform& t2) const;
        bool SwapsHandedness() const;
        explicit Transform(const Frame& frame);
        explicit Transform(Quaternion q);
        explicit operator Quaternion() const;

        void Decompose(Tangent3f* T, SquareMatrix<4>* R, SquareMatrix<4>* S) const;

        Point3fi operator()(const Point3fi& p) const {
            float x = float(p.x), y = float(p.y), z = float(p.z);
            // Compute transformed coordinates from point _x_, _y_, and _z_
            float xp = (m[0][0] * x + m[0][1] * y) + (m[0][2] * z + m[0][3]);
            float yp = (m[1][0] * x + m[1][1] * y) + (m[1][2] * z + m[1][3]);
            float zp = (m[2][0] * x + m[2][1] * y) + (m[2][2] * z + m[2][3]);
            float wp = (m[3][0] * x + m[3][1] * y) + (m[3][2] * z + m[3][3]);

            // Compute absolute error for transformed point, _pError_
            Tangent3f pError;
            if (p.IsExact()) {
                // Compute error for transformed exact _p_
                pError.x = gamma(3) * (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
                    std::abs(m[0][2] * z) + std::abs(m[0][3]));
                pError.y = gamma(3) * (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
                    std::abs(m[1][2] * z) + std::abs(m[1][3]));
                pError.z = gamma(3) * (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
                    std::abs(m[2][2] * z) + std::abs(m[2][3]));
            }
            else {
                // Compute error for transformed approximate _p_

                Tangent3f pInError = p.Error();
                pError.x = (gamma(3) + 1) * (std::abs(m[0][0]) * pInError.x +
                    std::abs(m[0][1]) * pInError.y +
                    std::abs(m[0][2]) * pInError.z) +
                    gamma(3) * (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
                        std::abs(m[0][2] * z) + std::abs(m[0][3]));
                pError.y = (gamma(3) + 1) * (std::abs(m[1][0]) * pInError.x +
                    std::abs(m[1][1]) * pInError.y +
                    std::abs(m[1][2]) * pInError.z) +
                    gamma(3) * (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
                        std::abs(m[1][2] * z) + std::abs(m[1][3]));
                pError.z = (gamma(3) + 1) * (std::abs(m[2][0]) * pInError.x +
                    std::abs(m[2][1]) * pInError.y +
                    std::abs(m[2][2]) * pInError.z) +
                    gamma(3) * (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
                        std::abs(m[2][2] * z) + std::abs(m[2][3]));
            }

            if (wp == 1)
                return Point3fi(Point3f(xp, yp, zp), pError);
            else
                return Point3fi(Point3f(xp / wp, yp / wp, zp / wp), pError / wp);
        }

        Tangent3fi operator()(const Tangent3fi& v) const;
        Point3fi ApplyInverse(const Point3fi& p) const;

    private:
        // Transform Private Members
        SquareMatrix<4> m, mInv;
    };

    // Transform Function Declarations

    Transform Translate(Tangent3f delta);
    Transform Scale(float x, float y, float z);
    Transform RotateX(float theta);
    Transform RotateY(float theta);
    Transform RotateZ(float theta);
    Transform LookAt(Point3f pos, Point3f look, Tangent3f up);
    Transform Orthographic(float znear, float zfar);
    Transform Perspective(float fov, float znear, float zfar);
    bool SolveLinearSystem2x2(const float A[2][2], const float B[2], float* x0, float* x1);

    // Transform Inline Functions
    inline Transform Inverse(const Transform& t) {
        return Transform(t.GetInverseMatrix(), t.GetMatrix());
    }

    inline Transform Transpose(const Transform& t) {
        return { Transpose(t.GetMatrix()), Transpose(t.GetInverseMatrix()) };
    }

    inline Transform Rotate(float sinTheta, float cosTheta, Tangent3f axis) {
        Tangent3f a = Normalize(axis);
        SquareMatrix<4> m;
        // Compute rotation of first basis vector
        m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
        m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
        m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
        m[0][3] = 0;

        // Compute rotations of second and third basis vectors
        m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
        m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
        m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
        m[1][3] = 0;

        m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
        m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
        m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
        m[2][3] = 0;

        return { m, Transpose(m) };
    }

    inline Transform Rotate(float theta, Tangent3f axis) {
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        return Rotate(sinTheta, cosTheta, axis);
    }

    inline Transform RotateFromTo(Tangent3f from, Tangent3f to) {
        // Compute intermediate vector for vector reflection
        Tangent3f refl;
        if (std::abs(from.x) < 0.72f && std::abs(to.x) < 0.72f)
            refl = Tangent3f(1, 0, 0);
        else if (std::abs(from.y) < 0.72f && std::abs(to.y) < 0.72f)
            refl = Tangent3f(0, 1, 0);
        else
            refl = Tangent3f(0, 0, 1);

        // Initialize matrix _r_ for rotation
        Tangent3f u = refl - from, v = refl - to;
        SquareMatrix<4> r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                // Initialize matrix element _r[i][j]_
                r[i][j] = ((i == j) ? 1 : 0) - 2 / Dot(u, u) * u[i] * u[j] -
                2 / Dot(v, v) * v[i] * v[j] +
                4 * Dot(u, v) / (Dot(u, u) * Dot(v, v)) * v[i] * u[j];

        return { r, Transpose(r) };
    }

    inline Tangent3fi Transform::operator()(const Tangent3fi& v) const {
        float x = float(v.x), y = float(v.y), z = float(v.z);
        Tangent3f vOutError;
        if (v.IsExact()) {
            vOutError.x = gamma(3) * (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
                std::abs(m[0][2] * z));
            vOutError.y = gamma(3) * (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
                std::abs(m[1][2] * z));
            vOutError.z = gamma(3) * (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
                std::abs(m[2][2] * z));
        }
        else {
            Tangent3f vInError = v.Error();
            vOutError.x = (gamma(3) + 1) * (std::abs(m[0][0]) * vInError.x +
                std::abs(m[0][1]) * vInError.y +
                std::abs(m[0][2]) * vInError.z) +
                gamma(3) * (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
                    std::abs(m[0][2] * z));
            vOutError.y = (gamma(3) + 1) * (std::abs(m[1][0]) * vInError.x +
                std::abs(m[1][1]) * vInError.y +
                std::abs(m[1][2]) * vInError.z) +
                gamma(3) * (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
                    std::abs(m[1][2] * z));
            vOutError.z = (gamma(3) + 1) * (std::abs(m[2][0]) * vInError.x +
                std::abs(m[2][1]) * vInError.y +
                std::abs(m[2][2]) * vInError.z) +
                gamma(3) * (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
                    std::abs(m[2][2] * z));
        }

        float xp = m[0][0] * x + m[0][1] * y + m[0][2] * z;
        float yp = m[1][0] * x + m[1][1] * y + m[1][2] * z;
        float zp = m[2][0] * x + m[2][1] * y + m[2][2] * z;

        return Tangent3fi(Tangent3f(xp, yp, zp), vOutError);
    }

    // Transform Inline Methods
    template <typename T>
    inline Point3<T> Transform::operator()(Point3<T> p) const {
        T xp = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
        T yp = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
        T zp = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
        T wp = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];
        if (wp == 1)
            return { xp, yp, zp };
        else
            return Point3<T>(xp / wp, yp / wp, zp / wp);
    }

    template <typename T>
    inline Tangent3<T> Transform::operator()(Tangent3<T> v) const {
        return { m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z };
    }

    template <typename T>
    inline Normal3<T> Transform::operator()(Normal3<T> n) const {
        T x = n.x, y = n.y, z = n.z;
        return { mInv[0][0] * x + mInv[1][0] * y + mInv[2][0] * z,
            mInv[0][1] * x + mInv[1][1] * y + mInv[2][1] * z,
            mInv[0][2] * x + mInv[1][2] * y + mInv[2][2] * z };
    }

    inline Ray Transform::operator()(const Ray& r, float* tMax) const {
        Point3fi o = (*this)(Point3fi(r.o));
        Tangent3f d = (*this)(r.d);
        // Offset ray origin to edge of error bounds and compute _tMax_
        if (float lengthSquared = LengthSquared(d); lengthSquared > 0) {
            float dt = Dot(Abs(d), o.Error()) / lengthSquared;
            o += d * dt;
            if (tMax)
                *tMax -= dt;
        }

        return Ray(Point3f(o), d);
    }

    inline RayDifferential Transform::operator()(const RayDifferential& r,
        float* tMax) const {
        Ray tr = (*this)(Ray(r), tMax);
        RayDifferential ret(tr.o, tr.d);
        ret.hasDifferentials = r.hasDifferentials;
        ret.rxOrigin = (*this)(r.rxOrigin);
        ret.ryOrigin = (*this)(r.ryOrigin);
        ret.rxDirection = (*this)(r.rxDirection);
        ret.ryDirection = (*this)(r.ryDirection);
        return ret;
    }


    template <typename T>
    inline Point3<T> Transform::operator()(const Point3<T>& p,
        Tangent3<T>* pError) const {
        T x = p.x, y = p.y, z = p.z;
        // Compute transformed coordinates from point _pt_
        T xp = (m[0][0] * x + m[0][1] * y) + (m[0][2] * z + m[0][3]);
        T yp = (m[1][0] * x + m[1][1] * y) + (m[1][2] * z + m[1][3]);
        T zp = (m[2][0] * x + m[2][1] * y) + (m[2][2] * z + m[2][3]);
        T wp = (m[3][0] * x + m[3][1] * y) + (m[3][2] * z + m[3][3]);

        // Compute absolute error for transformed point
        T xAbsSum = (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
            std::abs(m[0][2] * z) + std::abs(m[0][3]));
        T yAbsSum = (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
            std::abs(m[1][2] * z) + std::abs(m[1][3]));
        T zAbsSum = (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
            std::abs(m[2][2] * z) + std::abs(m[2][3]));
        *pError = gamma(3) * Tangent3<T>(xAbsSum, yAbsSum, zAbsSum);
        if (wp == 1)
            return Point3<T>(xp, yp, zp);
        else
            return Point3<T>(xp / wp, yp / wp, zp / wp);
    }

    template <typename T>
    inline Point3<T> Transform::operator()(const Point3<T>& pt,
        const Tangent3<T>& ptError,
        Tangent3<T>* absError) const {
        T x = pt.x, y = pt.y, z = pt.z;
        T xp = (m[0][0] * x + m[0][1] * y) + (m[0][2] * z + m[0][3]);
        T yp = (m[1][0] * x + m[1][1] * y) + (m[1][2] * z + m[1][3]);
        T zp = (m[2][0] * x + m[2][1] * y) + (m[2][2] * z + m[2][3]);
        T wp = (m[3][0] * x + m[3][1] * y) + (m[3][2] * z + m[3][3]);
        absError->x =
            (gamma(3) + (T)1) *
            (std::abs(m[0][0]) * ptError.x + std::abs(m[0][1]) * ptError.y +
                std::abs(m[0][2]) * ptError.z) +
            gamma(3) * (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
                std::abs(m[0][2] * z) + std::abs(m[0][3]));
        absError->y =
            (gamma(3) + (T)1) *
            (std::abs(m[1][0]) * ptError.x + std::abs(m[1][1]) * ptError.y +
                std::abs(m[1][2]) * ptError.z) +
            gamma(3) * (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
                std::abs(m[1][2] * z) + std::abs(m[1][3]));
        absError->z =
            (gamma(3) + (T)1) *
            (std::abs(m[2][0]) * ptError.x + std::abs(m[2][1]) * ptError.y +
                std::abs(m[2][2]) * ptError.z) +
            gamma(3) * (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
                std::abs(m[2][2] * z) + std::abs(m[2][3]));
        if (wp == 1.)
            return Point3<T>(xp, yp, zp);
        else
            return Point3<T>(xp / wp, yp / wp, zp / wp);
    }

    template <typename T>
    inline Tangent3<T> Transform::operator()(const Tangent3<T>& v,
        Tangent3<T>* absError) const {
        T x = v.x, y = v.y, z = v.z;
        absError->x =
            gamma(3) * (std::abs(m[0][0] * v.x) + std::abs(m[0][1] * v.y) +
                std::abs(m[0][2] * v.z));
        absError->y =
            gamma(3) * (std::abs(m[1][0] * v.x) + std::abs(m[1][1] * v.y) +
                std::abs(m[1][2] * v.z));
        absError->z =
            gamma(3) * (std::abs(m[2][0] * v.x) + std::abs(m[2][1] * v.y) +
                std::abs(m[2][2] * v.z));
        return Tangent3<T>(m[0][0] * x + m[0][1] * y + m[0][2] * z,
            m[1][0] * x + m[1][1] * y + m[1][2] * z,
            m[2][0] * x + m[2][1] * y + m[2][2] * z);
    }

    template <typename T>
    inline Tangent3<T> Transform::operator()(const Tangent3<T>& v,
        const Tangent3<T>& vError,
        Tangent3<T>* absError) const {
        T x = v.x, y = v.y, z = v.z;
        absError->x =
            (gamma(3) + (T)1) *
            (std::abs(m[0][0]) * vError.x + std::abs(m[0][1]) * vError.y +
                std::abs(m[0][2]) * vError.z) +
            gamma(3) * (std::abs(m[0][0] * v.x) + std::abs(m[0][1] * v.y) +
                std::abs(m[0][2] * v.z));
        absError->y =
            (gamma(3) + (T)1) *
            (std::abs(m[1][0]) * vError.x + std::abs(m[1][1]) * vError.y +
                std::abs(m[1][2]) * vError.z) +
            gamma(3) * (std::abs(m[1][0] * v.x) + std::abs(m[1][1] * v.y) +
                std::abs(m[1][2] * v.z));
        absError->z =
            (gamma(3) + (T)1) *
            (std::abs(m[2][0]) * vError.x + std::abs(m[2][1]) * vError.y +
                std::abs(m[2][2]) * vError.z) +
            gamma(3) * (std::abs(m[2][0] * v.x) + std::abs(m[2][1] * v.y) +
                std::abs(m[2][2] * v.z));
        return Tangent3<T>(m[0][0] * x + m[0][1] * y + m[0][2] * z,
            m[1][0] * x + m[1][1] * y + m[1][2] * z,
            m[2][0] * x + m[2][1] * y + m[2][2] * z);
    }

    inline Ray Transform::operator()(const Ray& r, Tangent3f* oError,
        Tangent3f* dError) const {
        Point3f o = (*this)(r.o, oError);
        Tangent3f d = (*this)(r.d, dError);
        float tMax = r.tMax;
        float lengthSquared = LengthSquared(d);
        if (lengthSquared > 0) {
            float dt = Dot(Abs(d), *oError) / lengthSquared;
            o += d * dt;
            tMax -= dt;
        }
        return Ray(o, d, tMax, r.time, r.medium);
    }

    inline Ray Transform::operator()(const Ray& r, const Tangent3f& oErrorIn,
        const Tangent3f& dErrorIn, Tangent3f* oErrorOut,
        Tangent3f* dErrorOut) const {
        Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
        Tangent3f d = (*this)(r.d, dErrorIn, dErrorOut);
        float tMax = r.tMax;
        float lengthSquared = LengthSquared(d);
        if (lengthSquared > 0) {
            float dt = Dot(Abs(d), *oErrorOut) / lengthSquared;
            o += d * dt;
            tMax -= dt;
        }
        return Ray(o, d, tMax, r.time, r.medium);
    }

    inline Transform::Transform(const Frame& frame)
        : Transform(SquareMatrix<4>(frame.x.x, frame.x.y, frame.x.z, 0.f, frame.y.x, frame.y.y,
            frame.y.z, 0.f, frame.z.x, frame.z.y, frame.z.z, 0.f, 0.f, 0.f, 0.f, 1.f)) {}

    inline Transform::Transform(Quaternion q) {
        float xx = q.v.x * q.v.x, yy = q.v.y * q.v.y, zz = q.v.z * q.v.z;
        float xy = q.v.x * q.v.y, xz = q.v.x * q.v.z, yz = q.v.y * q.v.z;
        float wx = q.v.x * q.w, wy = q.v.y * q.w, wz = q.v.z * q.w;

        m[0][0] = 1 - 2 * (yy + zz);
        m[0][1] = 2 * (xy + wz);
        m[0][2] = 2 * (xz - wy);
        m[1][0] = 2 * (xy - wz);
        m[1][1] = 1 - 2 * (xx + zz);
        m[1][2] = 2 * (yz + wx);
        m[2][0] = 2 * (xz + wy);
        m[2][1] = 2 * (yz - wx);
        m[2][2] = 1 - 2 * (xx + yy);

        // Transpose since we are left-handed.  Ugh.
        mInv = Transpose(m);
    }

    template <typename T>
    inline Point3<T> Transform::ApplyInverse(Point3<T> p) const {
        T x = p.x, y = p.y, z = p.z;
        T xp = (mInv[0][0] * x + mInv[0][1] * y) + (mInv[0][2] * z + mInv[0][3]);
        T yp = (mInv[1][0] * x + mInv[1][1] * y) + (mInv[1][2] * z + mInv[1][3]);
        T zp = (mInv[2][0] * x + mInv[2][1] * y) + (mInv[2][2] * z + mInv[2][3]);
        T wp = (mInv[3][0] * x + mInv[3][1] * y) + (mInv[3][2] * z + mInv[3][3]);
        if (wp == 1)
            return { xp, yp, zp };
        else
            return Point3<T>(xp / wp, yp / wp, zp / wp);
    }

    template <typename T>
    inline Tangent3<T> Transform::ApplyInverse(Tangent3<T> v) const {
        T x = v.x, y = v.y, z = v.z;
        return { mInv[0][0] * x + mInv[0][1] * y + mInv[0][2] * z,
            mInv[1][0] * x + mInv[1][1] * y + mInv[1][2] * z,
            mInv[2][0] * x + mInv[2][1] * y + mInv[2][2] * z };
    }

    template <typename T>
    inline Normal3<T> Transform::ApplyInverse(Normal3<T> n) const {
        T x = n.x, y = n.y, z = n.z;
        return { m[0][0] * x + m[1][0] * y + m[2][0] * z,
            m[0][1] * x + m[1][1] * y + m[2][1] * z,
            m[0][2] * x + m[1][2] * y + m[2][2] * z };
    }

    inline Ray Transform::ApplyInverse(const Ray& r, float* tMax) const {
        Point3fi o = ApplyInverse(Point3fi(r.o));
        Tangent3f d = ApplyInverse(r.d);
        // Offset ray origin to edge of error bounds and compute _tMax_
        float lengthSquared = LengthSquared(d);
        if (lengthSquared > 0) {
            Tangent3f oError(o.x.Width() / 2, o.y.Width() / 2, o.z.Width() / 2);
            float dt = Dot(Abs(d), oError) / lengthSquared;
            o += d * dt;
            if (tMax)
                *tMax -= dt;
        }
        return Ray(Point3f(o), d);
    }

    inline RayDifferential Transform::ApplyInverse(const RayDifferential& r,
        float* tMax) const {
        Ray tr = ApplyInverse(Ray(r), tMax);
        RayDifferential ret(tr.o, tr.d);
        ret.hasDifferentials = r.hasDifferentials;
        ret.rxOrigin = ApplyInverse(r.rxOrigin);
        ret.ryOrigin = ApplyInverse(r.ryOrigin);
        ret.rxDirection = ApplyInverse(r.rxDirection);
        ret.ryDirection = ApplyInverse(r.ryDirection);
        return ret;
    }

    class AnimatedTransform {
    public:
        // AnimatedTransform Public Methods
        AnimatedTransform(const Transform& startTransform, float startTime,
            const Transform& endTransform, float endTime);
        void Interpolate(float time, Transform& t) const;
        Ray operator()(const Ray& r) const;
        RayDifferential operator()(const RayDifferential& r) const;
        Point3f operator()(float time, const Point3f& p) const;
        Tangent3f operator()(float time, const Tangent3f& v) const;
        bool HasScale() const {
            return startTransform.HasScale() || endTransform.HasScale();
        }
        Bounds3f MotionBounds(const Bounds3f& b) const;
        Bounds3f BoundPointMotion(const Point3f& p) const;

    private:
        // AnimatedTransform Private Data
        const Transform startTransform, endTransform;
        const float startTime, endTime;
        const bool actuallyAnimated;
        Tangent3f T[2];
        Quaternion R[2];
        SquareMatrix<4> S[2];
        bool hasRotation;
        struct DerivativeTerm {
            DerivativeTerm() {}
            DerivativeTerm(float c, float x, float y, float z)
                : kc(c), kx(x), ky(y), kz(z) {}
            float kc, kx, ky, kz;
            float Eval(const Point3f& p) const {
                return kc + kx * p.x + ky * p.y + kz * p.z;
            }
        };
        DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
        static void FindZeros(float c1, float c2, float c3, float c4, float c5,
            float theta, Interval tInterval, float* zeros, int* nZeros, int depth = 8);
    };

} // namespace ligtfold