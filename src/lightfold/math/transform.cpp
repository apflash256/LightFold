#include <math/transform.h>

#include <core/interaction.h>

#include <iostream>

namespace lightfold {

    bool SolveLinearSystem2x2(const float A[2][2], const float B[2], float* x0, float* x1) {
        float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (std::abs(det) < 1e-10f) return false;
        *x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
        *x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
        if (std::isnan(*x0) || std::isnan(*x1)) return false;
        return true;
    }

    // Transform Function Definitions
    Transform Translate(Tangent3f delta) {
        SquareMatrix<4> m(1.f, 0.f, 0.f, delta.x,
            0.f, 1.f, 0.f, delta.y,
            0.f, 0.f, 1.f, delta.z,
            0.f, 0.f, 0.f, 1.f);
        SquareMatrix<4> minv(1.f, 0.f, 0.f, -delta.x,
            0.f, 1.f, 0.f, -delta.y,
            0.f, 0.f, 1.f, -delta.z,
            0.f, 0.f, 0.f, 1.f);
        return Transform(m, minv);
    }
    Transform Scale(float x, float y, float z) {
        SquareMatrix<4> m(x, 0.f, 0.f, 0.f,
            0.f, y, 0.f, 0.f,
            0.f, 0.f, z, 0.f,
            0.f, 0.f, 0.f, 1.f);
        SquareMatrix<4> minv(1 / x, 0.f, 0.f, 0.f,
            0.f, 1 / y, 0.f, 0.f,
            0.f, 0.f, 1 / z, 0.f,
            0.f, 0.f, 0.f, 1.f);
        return Transform(m, minv);
    }
    Transform RotateX(float theta) {
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        SquareMatrix<4> m(1.f, 0.f, 0.f, 0.f,
            0.f, cosTheta, -sinTheta, 0.f,
            0.f, sinTheta, cosTheta, 0.f,
            0.f, 0.f, 0.f, 1.f);
        return Transform(m, Transpose(m));
    }
    Transform RotateY(float theta) {
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        SquareMatrix<4> m(cosTheta, 0.f, sinTheta, 0.f,
            0.f, 1.f, 0.f, 0.f,
            -sinTheta, 0.f, cosTheta, 0.f,
            0.f, 0.f, 0.f, 1.f);
        return Transform(m, Transpose(m));
    }
    Transform RotateZ(float theta) {
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        SquareMatrix<4> m(cosTheta, -sinTheta, 0.f, 0.f,
            sinTheta, cosTheta, 0.f, 0.f,
            0.f, 0.f, 1.f, 0.f,
            0.f, 0.f, 0.f, 1.f);
        return Transform(m, Transpose(m));
    }

    Transform LookAt(Point3f pos, Point3f look, Tangent3f up) {
        SquareMatrix<4> CameraToWorld;
        // Initialize fourth column of viewing matrix
        CameraToWorld[0][3] = pos.x;
        CameraToWorld[1][3] = pos.y;
        CameraToWorld[2][3] = pos.z;
        CameraToWorld[3][3] = 1.f;

        // Initialize first three columns of viewing matrix
        Tangent3f dir = Normalize(look - pos);
        if (Length(Cross(Normalize(up), dir)) == 0)
            std::cout << "LookAt: \"up\" vector (" << up.x << ", " << up.y << ", " << up.z
            << ") and viewing direction (" << dir.x << ", " << dir.y << ", " << dir.z
            << " passed to LookAt are pointing in the same direction." << std::endl;
        Tangent3f right = Normalize(Cross(dir, Normalize(up)));
        Tangent3f newUp = Cross(right, dir);
        CameraToWorld[0][0] = -right.x;
        CameraToWorld[1][0] = -right.y;
        CameraToWorld[2][0] = -right.z;
        CameraToWorld[3][0] = 0.f;
        CameraToWorld[0][1] = newUp.x;
        CameraToWorld[1][1] = newUp.y;
        CameraToWorld[2][1] = newUp.z;
        CameraToWorld[3][1] = 0.f;
        CameraToWorld[0][2] = dir.x;
        CameraToWorld[1][2] = dir.y;
        CameraToWorld[2][2] = dir.z;
        CameraToWorld[3][2] = 0.f;

        return Transform(InvertOrExit(CameraToWorld), CameraToWorld);
    }

    Transform Orthographic(float near, float far) {
        // Map near to 0, far to 1
        return Scale(-1.f, -1.f, 1 / (far - near)) * Translate(Tangent3f(0.f, 0.f, near));
    }

    Transform Perspective(float fov, float near, float far) {
        // Perform projective divide for perspective projection
        SquareMatrix<4> persp(-1.f, 0.f, 0.f, 0.f,
            0.f, -1.f, 0.f, 0.f,
            0.f, 0.f, far / (far - near), -far * near / (far - near),
            0.f, 0.f, 1.f, 0.f);
        // Scale canonical perspective view to specified field of view
        float invTanAng = 1 / std::tan(fov / 2);
        return Scale(invTanAng, invTanAng, 1) * Transform(persp);
    }

    // Transform Method Definitions
    Bounds3f Transform::operator()(const Bounds3f& b) const {
        Bounds3f bt;
        for (int i = 0; i < 8; ++i)
            bt = Union(bt, (*this)(b.Corner(i)));
        return bt;
    }

    SurfaceInteraction Transform::operator()(const SurfaceInteraction& si) const {
        SurfaceInteraction ret;
        // Transform _p_ and _pError_ in _SurfaceInteraction_
        ret.p = (*this)(si.p, si.pError, &ret.pError);

        // Transform remaining members of _SurfaceInteraction_
        const Transform& t = *this;
        ret.n = Normalize(t(si.n));
        ret.wo = Normalize(t(si.wo));
        ret.time = si.time;
        ret.mediumInterface = si.mediumInterface;
        ret.uv = si.uv;
        ret.shape = si.shape;
        ret.dpdu = t(si.dpdu);
        ret.dpdv = t(si.dpdv);
        ret.dndu = t(si.dndu);
        ret.dndv = t(si.dndv);
        ret.shading.n = Normalize(t(si.shading.n));
        ret.shading.dpdu = t(si.shading.dpdu);
        ret.shading.dpdv = t(si.shading.dpdv);
        ret.shading.dndu = t(si.shading.dndu);
        ret.shading.dndv = t(si.shading.dndv);
        ret.dudx = si.dudx;
        ret.dvdx = si.dvdx;
        ret.dudy = si.dudy;
        ret.dvdy = si.dvdy;
        ret.dpdx = t(si.dpdx);
        ret.dpdy = t(si.dpdy);
        ret.bsdf = si.bsdf;
        ret.bssrdf = si.bssrdf;
        ret.primitive = si.primitive;
        //    ret.n = Faceforward(ret.n, ret.shading.n);
        ret.shading.n = FaceForward(ret.shading.n, ret.n);
        return ret;
    }

    Transform Transform::operator*(const Transform& t2) const {
        return Transform(m * t2.m, t2.mInv * mInv);
    }

    bool Transform::SwapsHandedness() const {
        // clang-format off
        SquareMatrix<3> s(m[0][0], m[0][1], m[0][2],
            m[1][0], m[1][1], m[1][2],
            m[2][0], m[2][1], m[2][2]);
        // clang-format on
        return Determinant(s) < 0;
    }

    Transform::operator Quaternion() const {
        float trace = m[0][0] + m[1][1] + m[2][2];
        Quaternion quat;
        if (trace > 0.f) {
            // Compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            float s = std::sqrt(trace + 1.0f);
            quat.w = s / 2.0f;
            s = 0.5f / s;
            quat.v.x = (m[2][1] - m[1][2]) * s;
            quat.v.y = (m[0][2] - m[2][0]) * s;
            quat.v.z = (m[1][0] - m[0][1]) * s;
        }
        else {
            // Compute largest of $x$, $y$, or $z$, then remaining components
            const int nxt[3] = { 1, 2, 0 };
            float q[3];
            int i = 0;
            if (m[1][1] > m[0][0])
                i = 1;
            if (m[2][2] > m[i][i])
                i = 2;
            int j = nxt[i];
            int k = nxt[j];
            float s = SafeSqrt((m[i][i] - (m[j][j] + m[k][k])) + 1.0f);
            q[i] = s * 0.5f;
            if (s != 0.f)
                s = 0.5f / s;
            quat.w = (m[k][j] - m[j][k]) * s;
            q[j] = (m[j][i] + m[i][j]) * s;
            q[k] = (m[k][i] + m[i][k]) * s;
            quat.v.x = q[0];
            quat.v.y = q[1];
            quat.v.z = q[2];
        }
        return quat;
    }

    void Transform::Decompose(Tangent3f* T, SquareMatrix<4>* R, SquareMatrix<4>* S) const {
        // Extract translation _T_ from transformation matrix
        T->x = m[0][3];
        T->y = m[1][3];
        T->z = m[2][3];

        // Compute new transformation matrix _M_ without translation
        SquareMatrix<4> M = m;
        for (int i = 0; i < 3; ++i)
            M[i][3] = M[3][i] = 0.f;
        M[3][3] = 1.f;

        // Extract rotation _R_ from transformation matrix
        float norm;
        int count = 0;
        *R = M;
        do {
            // Compute next matrix _Rnext_ in series
            SquareMatrix<4> Rit = InvertOrExit(Transpose(*R));
            SquareMatrix<4> Rnext = (*R + Rit) / 2;

            // Compute norm of difference between _R_ and _Rnext_
            norm = 0;
            for (int i = 0; i < 3; ++i) {
                float n = std::abs((*R)[i][0] - Rnext[i][0]) +
                    std::abs((*R)[i][1] - Rnext[i][1]) +
                    std::abs((*R)[i][2] - Rnext[i][2]);
                norm = std::max(norm, n);
            }

            *R = Rnext;
        } while (++count < 100 && norm > .0001);

        // Compute scale _S_ using rotation and original matrix
        *S = InvertOrExit(*R) * M;
    }

    AnimatedTransform::AnimatedTransform(const Transform& startTransform,
        float startTime,
        const Transform& endTransform,
        float endTime)
        : startTransform(startTransform),
        endTransform(endTransform),
        startTime(startTime),
        endTime(endTime),
        actuallyAnimated(startTransform != endTransform) {
        if (!actuallyAnimated)
            return;
        SquareMatrix<4> Rm;
        startTransform.Decompose(&T[0], &Rm, &S[0]);
        R[0] = Quaternion(Transform(Rm));
        endTransform.Decompose(&T[1], &Rm, &S[1]);
        R[1] = Quaternion(Transform(Rm));
        // Flip _R[1]_ if needed to select shortest path
        if (Dot(R[0], R[1]) < 0) R[1] = -R[1];
        hasRotation = Dot(R[0], R[1]) < 0.9995f;
        // Compute terms of motion derivative function
        if (hasRotation) {
            float cosTheta = Dot(R[0], R[1]);
            float theta = std::acos(Clamp(cosTheta, -1, 1));
            Quaternion qperp = Normalize(R[1] - R[0] * cosTheta);

            float t0x = T[0].x;
            float t0y = T[0].y;
            float t0z = T[0].z;
            float t1x = T[1].x;
            float t1y = T[1].y;
            float t1z = T[1].z;
            float q0x = R[0].v.x;
            float q0y = R[0].v.y;
            float q0z = R[0].v.z;
            float q0w = R[0].w;
            float qperpx = qperp.v.x;
            float qperpy = qperp.v.y;
            float qperpz = qperp.v.z;
            float qperpw = qperp.w;
            float s000 = S[0][0][0];
            float s001 = S[0][0][1];
            float s002 = S[0][0][2];
            float s010 = S[0][1][0];
            float s011 = S[0][1][1];
            float s012 = S[0][1][2];
            float s020 = S[0][2][0];
            float s021 = S[0][2][1];
            float s022 = S[0][2][2];
            float s100 = S[1][0][0];
            float s101 = S[1][0][1];
            float s102 = S[1][0][2];
            float s110 = S[1][1][0];
            float s111 = S[1][1][1];
            float s112 = S[1][1][2];
            float s120 = S[1][2][0];
            float s121 = S[1][2][1];
            float s122 = S[1][2][2];

            c1[0] = DerivativeTerm(
                -t0x + t1x,
                (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                s000 +
                q0w * q0z * s010 - qperpx * qperpy * s010 +
                qperpw * qperpz * s010 - q0w * q0y * s020 -
                qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
                q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
                qperpz * qperpz * s100 - q0w * q0z * s110 +
                qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                q0w * q0y * s120 + qperpw * qperpy * s120 +
                qperpx * qperpz * s120 +
                q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
                (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                s001 +
                q0w * q0z * s011 - qperpx * qperpy * s011 +
                qperpw * qperpz * s011 - q0w * q0y * s021 -
                qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
                q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
                qperpz * qperpz * s101 - q0w * q0z * s111 +
                qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                q0w * q0y * s121 + qperpw * qperpy * s121 +
                qperpx * qperpz * s121 +
                q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
                (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                s002 +
                q0w * q0z * s012 - qperpx * qperpy * s012 +
                qperpw * qperpz * s012 - q0w * q0y * s022 -
                qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
                q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
                qperpz * qperpz * s102 - q0w * q0z * s112 +
                qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                q0w * q0y * s122 + qperpw * qperpy * s122 +
                qperpx * qperpz * s122 +
                q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122));

            c2[0] = DerivativeTerm(
                0.,
                -(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
                qperpx * qperpy * s010 - qperpw * qperpz * s010 +
                qperpw * qperpy * s020 + qperpx * qperpz * s020 +
                q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
                qperpy * qperpy * s100 + qperpz * qperpz * s100 -
                qperpx * qperpy * s110 + qperpw * qperpz * s110 -
                qperpw * qperpy * s120 - qperpx * qperpz * s120 +
                2 * q0x * qperpy * s010 * theta -
                2 * q0w * qperpz * s010 * theta +
                2 * q0w * qperpy * s020 * theta +
                2 * q0x * qperpz * s020 * theta +
                q0y *
                (q0x * (-s010 + s110) + q0w * (-s020 + s120) +
                    2 * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020) *
                    theta) +
                q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) -
                    2 * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020) *
                    theta),
                -(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
                qperpx * qperpy * s011 - qperpw * qperpz * s011 +
                qperpw * qperpy * s021 + qperpx * qperpz * s021 +
                q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
                qperpy * qperpy * s101 + qperpz * qperpz * s101 -
                qperpx * qperpy * s111 + qperpw * qperpz * s111 -
                qperpw * qperpy * s121 - qperpx * qperpz * s121 +
                2 * q0x * qperpy * s011 * theta -
                2 * q0w * qperpz * s011 * theta +
                2 * q0w * qperpy * s021 * theta +
                2 * q0x * qperpz * s021 * theta +
                q0y *
                (q0x * (-s011 + s111) + q0w * (-s021 + s121) +
                    2 * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021) *
                    theta) +
                q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) -
                    2 * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021) *
                    theta),
                -(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
                qperpx * qperpy * s012 - qperpw * qperpz * s012 +
                qperpw * qperpy * s022 + qperpx * qperpz * s022 +
                q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
                qperpy * qperpy * s102 + qperpz * qperpz * s102 -
                qperpx * qperpy * s112 + qperpw * qperpz * s112 -
                qperpw * qperpy * s122 - qperpx * qperpz * s122 +
                2 * q0x * qperpy * s012 * theta -
                2 * q0w * qperpz * s012 * theta +
                2 * q0w * qperpy * s022 * theta +
                2 * q0x * qperpz * s022 * theta +
                q0y *
                (q0x * (-s012 + s112) + q0w * (-s022 + s122) +
                    2 * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022) *
                    theta) +
                q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) -
                    2 * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022) *
                    theta));

            c3[0] = DerivativeTerm(
                0.,
                -2 * (q0x * qperpy * s010 - q0w * qperpz * s010 +
                    q0w * qperpy * s020 + q0x * qperpz * s020 -
                    q0x * qperpy * s110 + q0w * qperpz * s110 -
                    q0w * qperpy * s120 - q0x * qperpz * s120 +
                    q0y * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                        2 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                    q0z * (-2 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                        2 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
                theta,
                -2 * (q0x * qperpy * s011 - q0w * qperpz * s011 +
                    q0w * qperpy * s021 + q0x * qperpz * s021 -
                    q0x * qperpy * s111 + q0w * qperpz * s111 -
                    q0w * qperpy * s121 - q0x * qperpz * s121 +
                    q0y * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                        2 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                    q0z * (-2 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                        2 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
                theta,
                -2 * (q0x * qperpy * s012 - q0w * qperpz * s012 +
                    q0w * qperpy * s022 + q0x * qperpz * s022 -
                    q0x * qperpy * s112 + q0w * qperpz * s112 -
                    q0w * qperpy * s122 - q0x * qperpz * s122 +
                    q0y * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                        2 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                    q0z * (-2 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                        2 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
                theta);

            c4[0] = DerivativeTerm(
                0.,
                -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
                q0x * qperpz * s020 + q0x * qperpy * s110 -
                q0w * qperpz * s110 + q0w * qperpy * s120 +
                q0x * qperpz * s120 + 2 * q0y * q0y * s000 * theta +
                2 * q0z * q0z * s000 * theta -
                2 * qperpy * qperpy * s000 * theta -
                2 * qperpz * qperpz * s000 * theta +
                2 * qperpx * qperpy * s010 * theta -
                2 * qperpw * qperpz * s010 * theta +
                2 * qperpw * qperpy * s020 * theta +
                2 * qperpx * qperpz * s020 * theta +
                q0y * (-(qperpx * s010) - qperpw * s020 +
                    2 * qperpy * (s000 - s100) + qperpx * s110 +
                    qperpw * s120 - 2 * q0x * s010 * theta -
                    2 * q0w * s020 * theta) +
                q0z * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020 -
                    2 * qperpz * s100 - qperpw * s110 + qperpx * s120 +
                    2 * q0w * s010 * theta - 2 * q0x * s020 * theta),
                -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
                q0x * qperpz * s021 + q0x * qperpy * s111 -
                q0w * qperpz * s111 + q0w * qperpy * s121 +
                q0x * qperpz * s121 + 2 * q0y * q0y * s001 * theta +
                2 * q0z * q0z * s001 * theta -
                2 * qperpy * qperpy * s001 * theta -
                2 * qperpz * qperpz * s001 * theta +
                2 * qperpx * qperpy * s011 * theta -
                2 * qperpw * qperpz * s011 * theta +
                2 * qperpw * qperpy * s021 * theta +
                2 * qperpx * qperpz * s021 * theta +
                q0y * (-(qperpx * s011) - qperpw * s021 +
                    2 * qperpy * (s001 - s101) + qperpx * s111 +
                    qperpw * s121 - 2 * q0x * s011 * theta -
                    2 * q0w * s021 * theta) +
                q0z * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021 -
                    2 * qperpz * s101 - qperpw * s111 + qperpx * s121 +
                    2 * q0w * s011 * theta - 2 * q0x * s021 * theta),
                -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
                q0x * qperpz * s022 + q0x * qperpy * s112 -
                q0w * qperpz * s112 + q0w * qperpy * s122 +
                q0x * qperpz * s122 + 2 * q0y * q0y * s002 * theta +
                2 * q0z * q0z * s002 * theta -
                2 * qperpy * qperpy * s002 * theta -
                2 * qperpz * qperpz * s002 * theta +
                2 * qperpx * qperpy * s012 * theta -
                2 * qperpw * qperpz * s012 * theta +
                2 * qperpw * qperpy * s022 * theta +
                2 * qperpx * qperpz * s022 * theta +
                q0y * (-(qperpx * s012) - qperpw * s022 +
                    2 * qperpy * (s002 - s102) + qperpx * s112 +
                    qperpw * s122 - 2 * q0x * s012 * theta -
                    2 * q0w * s022 * theta) +
                q0z * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022 -
                    2 * qperpz * s102 - qperpw * s112 + qperpx * s122 +
                    2 * q0w * s012 * theta - 2 * q0x * s022 * theta));

            c5[0] = DerivativeTerm(
                0.,
                2 * (qperpy * qperpy * s000 + qperpz * qperpz * s000 -
                    qperpx * qperpy * s010 + qperpw * qperpz * s010 -
                    qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                    qperpy * qperpy * s100 - qperpz * qperpz * s100 +
                    q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) +
                    qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                    q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) +
                    qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                    q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
                theta,
                2 * (qperpy * qperpy * s001 + qperpz * qperpz * s001 -
                    qperpx * qperpy * s011 + qperpw * qperpz * s011 -
                    qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                    qperpy * qperpy * s101 - qperpz * qperpz * s101 +
                    q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) +
                    qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                    q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) +
                    qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                    q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
                theta,
                2 * (qperpy * qperpy * s002 + qperpz * qperpz * s002 -
                    qperpx * qperpy * s012 + qperpw * qperpz * s012 -
                    qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                    qperpy * qperpy * s102 - qperpz * qperpz * s102 +
                    q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) +
                    qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                    q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) +
                    qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                    q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
                theta);

            c1[1] = DerivativeTerm(
                -t0y + t1y,
                -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
                q0z * q0z * s010 + qperpx * qperpx * s010 +
                qperpz * qperpz * s010 - q0y * q0z * s020 +
                qperpw * qperpx * s020 - qperpy * qperpz * s020 +
                qperpx * qperpy * s100 + qperpw * qperpz * s100 +
                q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
                q0z * q0z * s110 - qperpx * qperpx * s110 -
                qperpz * qperpz * s110 +
                q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                q0y * q0z * s120 - qperpw * qperpx * s120 +
                qperpy * qperpz * s120,
                -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
                q0z * q0z * s011 + qperpx * qperpx * s011 +
                qperpz * qperpz * s011 - q0y * q0z * s021 +
                qperpw * qperpx * s021 - qperpy * qperpz * s021 +
                qperpx * qperpy * s101 + qperpw * qperpz * s101 +
                q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
                q0z * q0z * s111 - qperpx * qperpx * s111 -
                qperpz * qperpz * s111 +
                q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                q0y * q0z * s121 - qperpw * qperpx * s121 +
                qperpy * qperpz * s121,
                -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
                q0z * q0z * s012 + qperpx * qperpx * s012 +
                qperpz * qperpz * s012 - q0y * q0z * s022 +
                qperpw * qperpx * s022 - qperpy * qperpz * s022 +
                qperpx * qperpy * s102 + qperpw * qperpz * s102 +
                q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
                q0z * q0z * s112 - qperpx * qperpx * s112 -
                qperpz * qperpz * s112 +
                q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                q0y * q0z * s122 - qperpw * qperpx * s122 +
                qperpy * qperpz * s122);

            c2[1] = DerivativeTerm(
                0.,
                qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
                qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                q0y * q0z * s020 - qperpw * qperpx * s020 +
                qperpy * qperpz * s020 - qperpx * qperpy * s100 -
                qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
                q0z * q0z * s110 + qperpx * qperpx * s110 +
                qperpz * qperpz * s110 + q0y * q0z * s120 +
                qperpw * qperpx * s120 - qperpy * qperpz * s120 +
                2 * q0z * qperpw * s000 * theta +
                2 * q0y * qperpx * s000 * theta -
                4 * q0z * qperpz * s010 * theta +
                2 * q0z * qperpy * s020 * theta +
                2 * q0y * qperpz * s020 * theta +
                q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 +
                    2 * qperpy * s000 * theta - 4 * qperpx * s010 * theta -
                    2 * qperpw * s020 * theta) +
                q0w * (-(q0z * s000) + q0z * s100 + 2 * qperpz * s000 * theta -
                    2 * qperpx * s020 * theta),
                qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
                qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                q0y * q0z * s021 - qperpw * qperpx * s021 +
                qperpy * qperpz * s021 - qperpx * qperpy * s101 -
                qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
                q0z * q0z * s111 + qperpx * qperpx * s111 +
                qperpz * qperpz * s111 + q0y * q0z * s121 +
                qperpw * qperpx * s121 - qperpy * qperpz * s121 +
                2 * q0z * qperpw * s001 * theta +
                2 * q0y * qperpx * s001 * theta -
                4 * q0z * qperpz * s011 * theta +
                2 * q0z * qperpy * s021 * theta +
                2 * q0y * qperpz * s021 * theta +
                q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 +
                    2 * qperpy * s001 * theta - 4 * qperpx * s011 * theta -
                    2 * qperpw * s021 * theta) +
                q0w * (-(q0z * s001) + q0z * s101 + 2 * qperpz * s001 * theta -
                    2 * qperpx * s021 * theta),
                qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
                qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                q0y * q0z * s022 - qperpw * qperpx * s022 +
                qperpy * qperpz * s022 - qperpx * qperpy * s102 -
                qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
                q0z * q0z * s112 + qperpx * qperpx * s112 +
                qperpz * qperpz * s112 + q0y * q0z * s122 +
                qperpw * qperpx * s122 - qperpy * qperpz * s122 +
                2 * q0z * qperpw * s002 * theta +
                2 * q0y * qperpx * s002 * theta -
                4 * q0z * qperpz * s012 * theta +
                2 * q0z * qperpy * s022 * theta +
                2 * q0y * qperpz * s022 * theta +
                q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 +
                    2 * qperpy * s002 * theta - 4 * qperpx * s012 * theta -
                    2 * qperpw * s022 * theta) +
                q0w * (-(q0z * s002) + q0z * s102 + 2 * qperpz * s002 * theta -
                    2 * qperpx * s022 * theta));

            c3[1] = DerivativeTerm(
                0., 2 * (-(q0x * qperpy * s000) - q0w * qperpz * s000 +
                    2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
                    q0w * qperpx * s020 + q0x * qperpy * s100 +
                    q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
                    q0x * qperpw * s120 - q0w * qperpx * s120 +
                    q0z * (2 * qperpz * s010 - qperpy * s020 +
                        qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                        qperpy * s120) +
                    q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                        qperpz * s120)) *
                theta,
                2 * (-(q0x * qperpy * s001) - q0w * qperpz * s001 +
                    2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
                    q0w * qperpx * s021 + q0x * qperpy * s101 +
                    q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
                    q0x * qperpw * s121 - q0w * qperpx * s121 +
                    q0z * (2 * qperpz * s011 - qperpy * s021 +
                        qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                        qperpy * s121) +
                    q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                        qperpz * s121)) *
                theta,
                2 * (-(q0x * qperpy * s002) - q0w * qperpz * s002 +
                    2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
                    q0w * qperpx * s022 + q0x * qperpy * s102 +
                    q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
                    q0x * qperpw * s122 - q0w * qperpx * s122 +
                    q0z * (2 * qperpz * s012 - qperpy * s022 +
                        qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                        qperpy * s122) +
                    q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                        qperpz * s122)) *
                theta);

            c4[1] = DerivativeTerm(
                0.,
                -(q0x * qperpy * s000) - q0w * qperpz * s000 +
                2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
                q0w * qperpx * s020 + q0x * qperpy * s100 +
                q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
                q0x * qperpw * s120 - q0w * qperpx * s120 +
                2 * qperpx * qperpy * s000 * theta +
                2 * qperpw * qperpz * s000 * theta +
                2 * q0x * q0x * s010 * theta + 2 * q0z * q0z * s010 * theta -
                2 * qperpx * qperpx * s010 * theta -
                2 * qperpz * qperpz * s010 * theta +
                2 * q0w * q0x * s020 * theta -
                2 * qperpw * qperpx * s020 * theta +
                2 * qperpy * qperpz * s020 * theta +
                q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                    qperpz * s120 - 2 * q0x * s000 * theta) +
                q0z * (2 * qperpz * s010 - qperpy * s020 +
                    qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                    qperpy * s120 - 2 * q0w * s000 * theta -
                    2 * q0y * s020 * theta),
                -(q0x * qperpy * s001) - q0w * qperpz * s001 +
                2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
                q0w * qperpx * s021 + q0x * qperpy * s101 +
                q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
                q0x * qperpw * s121 - q0w * qperpx * s121 +
                2 * qperpx * qperpy * s001 * theta +
                2 * qperpw * qperpz * s001 * theta +
                2 * q0x * q0x * s011 * theta + 2 * q0z * q0z * s011 * theta -
                2 * qperpx * qperpx * s011 * theta -
                2 * qperpz * qperpz * s011 * theta +
                2 * q0w * q0x * s021 * theta -
                2 * qperpw * qperpx * s021 * theta +
                2 * qperpy * qperpz * s021 * theta +
                q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                    qperpz * s121 - 2 * q0x * s001 * theta) +
                q0z * (2 * qperpz * s011 - qperpy * s021 +
                    qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                    qperpy * s121 - 2 * q0w * s001 * theta -
                    2 * q0y * s021 * theta),
                -(q0x * qperpy * s002) - q0w * qperpz * s002 +
                2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
                q0w * qperpx * s022 + q0x * qperpy * s102 +
                q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
                q0x * qperpw * s122 - q0w * qperpx * s122 +
                2 * qperpx * qperpy * s002 * theta +
                2 * qperpw * qperpz * s002 * theta +
                2 * q0x * q0x * s012 * theta + 2 * q0z * q0z * s012 * theta -
                2 * qperpx * qperpx * s012 * theta -
                2 * qperpz * qperpz * s012 * theta +
                2 * q0w * q0x * s022 * theta -
                2 * qperpw * qperpx * s022 * theta +
                2 * qperpy * qperpz * s022 * theta +
                q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                    qperpz * s122 - 2 * q0x * s002 * theta) +
                q0z * (2 * qperpz * s012 - qperpy * s022 +
                    qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                    qperpy * s122 - 2 * q0w * s002 * theta -
                    2 * q0y * s022 * theta));

            c5[1] = DerivativeTerm(
                0., -2 * (qperpx * qperpy * s000 + qperpw * qperpz * s000 +
                    q0z * q0z * s010 - qperpx * qperpx * s010 -
                    qperpz * qperpz * s010 - q0y * q0z * s020 -
                    qperpw * qperpx * s020 + qperpy * qperpz * s020 -
                    qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                    q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) -
                    q0z * q0z * s110 + qperpx * qperpx * s110 +
                    qperpz * qperpz * s110 +
                    q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                    q0y * q0z * s120 + qperpw * qperpx * s120 -
                    qperpy * qperpz * s120) *
                theta,
                -2 * (qperpx * qperpy * s001 + qperpw * qperpz * s001 +
                    q0z * q0z * s011 - qperpx * qperpx * s011 -
                    qperpz * qperpz * s011 - q0y * q0z * s021 -
                    qperpw * qperpx * s021 + qperpy * qperpz * s021 -
                    qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                    q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) -
                    q0z * q0z * s111 + qperpx * qperpx * s111 +
                    qperpz * qperpz * s111 +
                    q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                    q0y * q0z * s121 + qperpw * qperpx * s121 -
                    qperpy * qperpz * s121) *
                theta,
                -2 * (qperpx * qperpy * s002 + qperpw * qperpz * s002 +
                    q0z * q0z * s012 - qperpx * qperpx * s012 -
                    qperpz * qperpz * s012 - q0y * q0z * s022 -
                    qperpw * qperpx * s022 + qperpy * qperpz * s022 -
                    qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                    q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) -
                    q0z * q0z * s112 + qperpx * qperpx * s112 +
                    qperpz * qperpz * s112 +
                    q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                    q0y * q0z * s122 + qperpw * qperpx * s122 -
                    qperpy * qperpz * s122) *
                theta);

            c1[2] = DerivativeTerm(
                -t0z + t1z, (qperpw * qperpy * s000 - qperpx * qperpz * s000 -
                    q0y * q0z * s010 - qperpw * qperpx * s010 -
                    qperpy * qperpz * s010 - s020 + q0y * q0y * s020 +
                    qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                    qperpw * qperpy * s100 + qperpx * qperpz * s100 +
                    q0x * q0z * (-s000 + s100) + q0y * q0z * s110 +
                    qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                    q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) +
                    q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 -
                    qperpx * qperpx * s120 - qperpy * qperpy * s120),
                (qperpw * qperpy * s001 - qperpx * qperpz * s001 -
                    q0y * q0z * s011 - qperpw * qperpx * s011 -
                    qperpy * qperpz * s011 - s021 + q0y * q0y * s021 +
                    qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                    qperpw * qperpy * s101 + qperpx * qperpz * s101 +
                    q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
                    qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                    q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
                    q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
                    qperpx * qperpx * s121 - qperpy * qperpy * s121),
                (qperpw * qperpy * s002 - qperpx * qperpz * s002 -
                    q0y * q0z * s012 - qperpw * qperpx * s012 -
                    qperpy * qperpz * s012 - s022 + q0y * q0y * s022 +
                    qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                    qperpw * qperpy * s102 + qperpx * qperpz * s102 +
                    q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
                    qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                    q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
                    q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
                    qperpx * qperpx * s122 - qperpy * qperpy * s122));

            c2[2] = DerivativeTerm(
                0.,
                (q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
                    qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
                    qperpw * qperpx * s010 + qperpy * qperpz * s010 +
                    q0x * q0x * s020 + q0y * q0y * s020 - qperpx * qperpx * s020 -
                    qperpy * qperpy * s020 - q0w * q0y * s100 + q0x * q0z * s100 +
                    qperpw * qperpy * s100 - qperpx * qperpz * s100 +
                    q0w * q0x * s110 + q0y * q0z * s110 - qperpw * qperpx * s110 -
                    qperpy * qperpz * s110 - q0x * q0x * s120 - q0y * q0y * s120 +
                    qperpx * qperpx * s120 + qperpy * qperpy * s120 -
                    2 * q0y * qperpw * s000 * theta + 2 * q0z * qperpx * s000 * theta -
                    2 * q0w * qperpy * s000 * theta + 2 * q0x * qperpz * s000 * theta +
                    2 * q0x * qperpw * s010 * theta + 2 * q0w * qperpx * s010 * theta +
                    2 * q0z * qperpy * s010 * theta + 2 * q0y * qperpz * s010 * theta -
                    4 * q0x * qperpx * s020 * theta - 4 * q0y * qperpy * s020 * theta),
                (q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
                    qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
                    qperpw * qperpx * s011 + qperpy * qperpz * s011 +
                    q0x * q0x * s021 + q0y * q0y * s021 - qperpx * qperpx * s021 -
                    qperpy * qperpy * s021 - q0w * q0y * s101 + q0x * q0z * s101 +
                    qperpw * qperpy * s101 - qperpx * qperpz * s101 +
                    q0w * q0x * s111 + q0y * q0z * s111 - qperpw * qperpx * s111 -
                    qperpy * qperpz * s111 - q0x * q0x * s121 - q0y * q0y * s121 +
                    qperpx * qperpx * s121 + qperpy * qperpy * s121 -
                    2 * q0y * qperpw * s001 * theta + 2 * q0z * qperpx * s001 * theta -
                    2 * q0w * qperpy * s001 * theta + 2 * q0x * qperpz * s001 * theta +
                    2 * q0x * qperpw * s011 * theta + 2 * q0w * qperpx * s011 * theta +
                    2 * q0z * qperpy * s011 * theta + 2 * q0y * qperpz * s011 * theta -
                    4 * q0x * qperpx * s021 * theta - 4 * q0y * qperpy * s021 * theta),
                (q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
                    qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
                    qperpw * qperpx * s012 + qperpy * qperpz * s012 +
                    q0x * q0x * s022 + q0y * q0y * s022 - qperpx * qperpx * s022 -
                    qperpy * qperpy * s022 - q0w * q0y * s102 + q0x * q0z * s102 +
                    qperpw * qperpy * s102 - qperpx * qperpz * s102 +
                    q0w * q0x * s112 + q0y * q0z * s112 - qperpw * qperpx * s112 -
                    qperpy * qperpz * s112 - q0x * q0x * s122 - q0y * q0y * s122 +
                    qperpx * qperpx * s122 + qperpy * qperpy * s122 -
                    2 * q0y * qperpw * s002 * theta + 2 * q0z * qperpx * s002 * theta -
                    2 * q0w * qperpy * s002 * theta + 2 * q0x * qperpz * s002 * theta +
                    2 * q0x * qperpw * s012 * theta + 2 * q0w * qperpx * s012 * theta +
                    2 * q0z * qperpy * s012 * theta + 2 * q0y * qperpz * s012 * theta -
                    4 * q0x * qperpx * s022 * theta -
                    4 * q0y * qperpy * s022 * theta));

            c3[2] = DerivativeTerm(
                0., -2 * (-(q0w * qperpy * s000) + q0x * qperpz * s000 +
                    q0x * qperpw * s010 + q0w * qperpx * s010 -
                    2 * q0x * qperpx * s020 + q0w * qperpy * s100 -
                    q0x * qperpz * s100 - q0x * qperpw * s110 -
                    q0w * qperpx * s110 +
                    q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 -
                        qperpy * s110) +
                    2 * q0x * qperpx * s120 +
                    q0y * (qperpz * s010 - 2 * qperpy * s020 +
                        qperpw * (-s000 + s100) - qperpz * s110 +
                        2 * qperpy * s120)) *
                theta,
                -2 * (-(q0w * qperpy * s001) + q0x * qperpz * s001 +
                    q0x * qperpw * s011 + q0w * qperpx * s011 -
                    2 * q0x * qperpx * s021 + q0w * qperpy * s101 -
                    q0x * qperpz * s101 - q0x * qperpw * s111 -
                    q0w * qperpx * s111 +
                    q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 -
                        qperpy * s111) +
                    2 * q0x * qperpx * s121 +
                    q0y * (qperpz * s011 - 2 * qperpy * s021 +
                        qperpw * (-s001 + s101) - qperpz * s111 +
                        2 * qperpy * s121)) *
                theta,
                -2 * (-(q0w * qperpy * s002) + q0x * qperpz * s002 +
                    q0x * qperpw * s012 + q0w * qperpx * s012 -
                    2 * q0x * qperpx * s022 + q0w * qperpy * s102 -
                    q0x * qperpz * s102 - q0x * qperpw * s112 -
                    q0w * qperpx * s112 +
                    q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 -
                        qperpy * s112) +
                    2 * q0x * qperpx * s122 +
                    q0y * (qperpz * s012 - 2 * qperpy * s022 +
                        qperpw * (-s002 + s102) - qperpz * s112 +
                        2 * qperpy * s122)) *
                theta);

            c4[2] = DerivativeTerm(
                0.,
                q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
                q0w * qperpx * s010 + 2 * q0x * qperpx * s020 -
                q0w * qperpy * s100 + q0x * qperpz * s100 +
                q0x * qperpw * s110 + q0w * qperpx * s110 -
                2 * q0x * qperpx * s120 - 2 * qperpw * qperpy * s000 * theta +
                2 * qperpx * qperpz * s000 * theta -
                2 * q0w * q0x * s010 * theta +
                2 * qperpw * qperpx * s010 * theta +
                2 * qperpy * qperpz * s010 * theta +
                2 * q0x * q0x * s020 * theta + 2 * q0y * q0y * s020 * theta -
                2 * qperpx * qperpx * s020 * theta -
                2 * qperpy * qperpy * s020 * theta +
                q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 +
                    qperpy * s110 - 2 * q0x * s000 * theta) +
                q0y * (-(qperpz * s010) + 2 * qperpy * s020 +
                    qperpw * (s000 - s100) + qperpz * s110 -
                    2 * qperpy * s120 + 2 * q0w * s000 * theta -
                    2 * q0z * s010 * theta),
                q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
                q0w * qperpx * s011 + 2 * q0x * qperpx * s021 -
                q0w * qperpy * s101 + q0x * qperpz * s101 +
                q0x * qperpw * s111 + q0w * qperpx * s111 -
                2 * q0x * qperpx * s121 - 2 * qperpw * qperpy * s001 * theta +
                2 * qperpx * qperpz * s001 * theta -
                2 * q0w * q0x * s011 * theta +
                2 * qperpw * qperpx * s011 * theta +
                2 * qperpy * qperpz * s011 * theta +
                2 * q0x * q0x * s021 * theta + 2 * q0y * q0y * s021 * theta -
                2 * qperpx * qperpx * s021 * theta -
                2 * qperpy * qperpy * s021 * theta +
                q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 +
                    qperpy * s111 - 2 * q0x * s001 * theta) +
                q0y * (-(qperpz * s011) + 2 * qperpy * s021 +
                    qperpw * (s001 - s101) + qperpz * s111 -
                    2 * qperpy * s121 + 2 * q0w * s001 * theta -
                    2 * q0z * s011 * theta),
                q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
                q0w * qperpx * s012 + 2 * q0x * qperpx * s022 -
                q0w * qperpy * s102 + q0x * qperpz * s102 +
                q0x * qperpw * s112 + q0w * qperpx * s112 -
                2 * q0x * qperpx * s122 - 2 * qperpw * qperpy * s002 * theta +
                2 * qperpx * qperpz * s002 * theta -
                2 * q0w * q0x * s012 * theta +
                2 * qperpw * qperpx * s012 * theta +
                2 * qperpy * qperpz * s012 * theta +
                2 * q0x * q0x * s022 * theta + 2 * q0y * q0y * s022 * theta -
                2 * qperpx * qperpx * s022 * theta -
                2 * qperpy * qperpy * s022 * theta +
                q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 +
                    qperpy * s112 - 2 * q0x * s002 * theta) +
                q0y * (-(qperpz * s012) + 2 * qperpy * s022 +
                    qperpw * (s002 - s102) + qperpz * s112 -
                    2 * qperpy * s122 + 2 * q0w * s002 * theta -
                    2 * q0z * s012 * theta));

            c5[2] = DerivativeTerm(
                0., 2 * (qperpw * qperpy * s000 - qperpx * qperpz * s000 +
                    q0y * q0z * s010 - qperpw * qperpx * s010 -
                    qperpy * qperpz * s010 - q0y * q0y * s020 +
                    qperpx * qperpx * s020 + qperpy * qperpy * s020 +
                    q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 +
                    qperpx * qperpz * s100 +
                    q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) -
                    q0y * q0z * s110 + qperpw * qperpx * s110 +
                    qperpy * qperpz * s110 + q0y * q0y * s120 -
                    qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                    q0x * q0x * (-s020 + s120)) *
                theta,
                2 * (qperpw * qperpy * s001 - qperpx * qperpz * s001 +
                    q0y * q0z * s011 - qperpw * qperpx * s011 -
                    qperpy * qperpz * s011 - q0y * q0y * s021 +
                    qperpx * qperpx * s021 + qperpy * qperpy * s021 +
                    q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 +
                    qperpx * qperpz * s101 +
                    q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) -
                    q0y * q0z * s111 + qperpw * qperpx * s111 +
                    qperpy * qperpz * s111 + q0y * q0y * s121 -
                    qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                    q0x * q0x * (-s021 + s121)) *
                theta,
                2 * (qperpw * qperpy * s002 - qperpx * qperpz * s002 +
                    q0y * q0z * s012 - qperpw * qperpx * s012 -
                    qperpy * qperpz * s012 - q0y * q0y * s022 +
                    qperpx * qperpx * s022 + qperpy * qperpy * s022 +
                    q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 +
                    qperpx * qperpz * s102 +
                    q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) -
                    q0y * q0z * s112 + qperpw * qperpx * s112 +
                    qperpy * qperpz * s112 + q0y * q0y * s122 -
                    qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                    q0x * q0x * (-s022 + s122)) *
                theta);
        }
    }

    void AnimatedTransform::Interpolate(float time, Transform &t) const {
        // Handle boundary conditions for matrix interpolation
        if (!actuallyAnimated || time <= startTime) {
            t = startTransform;
            return;
        }
        if (time >= endTime) {
            t = endTransform;
            return;
        }
        float dt = (time - startTime) / (endTime - startTime);
        // Interpolate translation at _dt_
        Tangent3f trans = (1 - dt) * T[0] + dt * T[1];

        // Interpolate rotation at _dt_
        Quaternion rotate = Slerp(dt, R[0], R[1]);

        // Interpolate scale at _dt_
        SquareMatrix<4> scale = (1 - dt) * S[0] + dt * S[1];

        // Compute interpolated matrix as product of interpolated components
        t = Translate(trans) * Transform(rotate) * Transform(scale);
    }

    Ray AnimatedTransform::operator()(const Ray& r) const {
        if (!actuallyAnimated || r.time <= startTime)
            return (startTransform)(r);
        else if (r.time >= endTime)
            return (endTransform)(r);
        else {
            Transform t;
            Interpolate(r.time, t);
            return t(r);
        }
    }

    RayDifferential AnimatedTransform::operator()(const RayDifferential& r) const {
        if (!actuallyAnimated || r.time <= startTime)
            return (startTransform)(r);
        else if (r.time >= endTime)
            return (endTransform)(r);
        else {
            Transform t;
            Interpolate(r.time, t);
            return t(r);
        }
    }

    Point3f AnimatedTransform::operator()(float time, const Point3f& p) const {
        if (!actuallyAnimated || time <= startTime)
            return (startTransform)(p);
        else if (time >= endTime)
            return (endTransform)(p);
        Transform t;
        Interpolate(time, t);
        return t(p);
    }

    Tangent3f AnimatedTransform::operator()(float time, const Tangent3f& v) const {
        if (!actuallyAnimated || time <= startTime)
            return (startTransform)(v);
        else if (time >= endTime)
            return (endTransform)(v);
        Transform t;
        Interpolate(time, t);
        return t(v);
    }

    Bounds3f AnimatedTransform::MotionBounds(const Bounds3f& b) const {
        if (!actuallyAnimated) return (startTransform)(b);
        if (hasRotation == false) return Union((startTransform)(b), (endTransform)(b));
        // Return motion bounds accounting for animated rotation
        Bounds3f bounds;
        for (int corner = 0; corner < 8; ++corner)
            bounds = Union(bounds, BoundPointMotion(b.Corner(corner)));
        return bounds;
    }

    Bounds3f AnimatedTransform::BoundPointMotion(const Point3f& p) const {
        if (!actuallyAnimated) return Bounds3f((startTransform)(p));
        Bounds3f bounds((startTransform)(p), (endTransform)(p));
        float cosTheta = Dot(R[0], R[1]);
        float theta = std::acos(Clamp(cosTheta, -1, 1));
        for (int c = 0; c < 3; ++c) {
            // Find any motion derivative zeros for the component _c_
            float zeros[8];
            int nZeros = 0;
            FindZeros(c1[c].Eval(p), c2[c].Eval(p), c3[c].Eval(p),
                c4[c].Eval(p), c5[c].Eval(p), theta, Interval(0., 1.),
                zeros, &nZeros);
            // Expand bounding box for any motion derivative zeros found
            for (int i = 0; i < nZeros; ++i) {
                Point3f pz = (*this)(std::lerp(zeros[i], startTime, endTime), p);
                bounds = Union(bounds, pz);
            }
        }
        return bounds;
    }

    void AnimatedTransform::FindZeros(float c1, float c2, float c3, float c4, float c5,
        float theta, Interval tInterval,
        float* zeros, int* nZeros, int depth) {
        // Evaluate motion derivative in interval form, return if no zeros
        Interval dadt =
            Interval(c1) +
            (Interval(c2) + Interval(c3) * tInterval) * Cos(Interval(2 * theta) * tInterval) +
            (Interval(c4) + Interval(c5) * tInterval) * Sin(Interval(2 * theta) * tInterval);
        if (dadt.LowerBound() > 0 || dadt.UpperBound() < 0 ||
            dadt.LowerBound() == dadt.UpperBound())
            return;

        // Either split range and recurse or report a zero
        if (depth > 0 && dadt.Width() > 1e-3) {
            // Split _tInterval_ and check both resulting intervals
            float mid = tInterval.Midpoint();
            FindZeros(c1, c2, c3, c4, c5, theta, Interval(tInterval.LowerBound(), mid), zeros,
                nZeros, depth - 1);
            FindZeros(c1, c2, c3, c4, c5, theta, Interval(mid, tInterval.UpperBound()), zeros,
                nZeros, depth - 1);

        }
        else {
            // Use Newton's method to refine zero
            float tNewton = tInterval.Midpoint();
            for (int i = 0; i < 4; ++i) {
                // Evaluate motion function derivative and its derivative at _tNewton_
                float fNewton = c1 + (c2 + c3 * tNewton) * std::cos(2 * theta * tNewton) +
                    (c4 + c5 * tNewton) * std::sin(2 * theta * tNewton);
                float fPrimeNewton =
                    (c3 + 2 * (c4 + c5 * tNewton) * theta) * std::cos(2 * tNewton * theta) +
                    (c5 - 2 * (c2 + c3 * tNewton) * theta) * std::sin(2 * tNewton * theta);

                if (fNewton == 0 || fPrimeNewton == 0)
                    break;
                tNewton = tNewton - fNewton / fPrimeNewton;
            }

            // Record zero if refined zero is in interval
            if (tNewton >= tInterval.LowerBound() - 1e-3f &&
                tNewton < tInterval.UpperBound() + 1e-3f) {
                zeros[*nZeros] = tNewton;
                (*nZeros)++;
            }
        }
    }

} // namespace lightfold
