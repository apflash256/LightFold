#include <math/transform.h>

#include <iostream>

namespace lightfold {

    // Transform Function Definitions
    Transform Translate(TanVector3f delta) {
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

    Transform LookAt(Point3f pos, Point3f look, TanVector3f up) {
        SquareMatrix<4> worldFromCamera;
        // Initialize fourth column of viewing matrix
        worldFromCamera[0][3] = pos.x;
        worldFromCamera[1][3] = pos.y;
        worldFromCamera[2][3] = pos.z;
        worldFromCamera[3][3] = 1.f;

        // Initialize first three columns of viewing matrix
        TanVector3f dir = Normalize(look - pos);
        if (Length(Cross(Normalize(up), dir)) == 0)
            std::cout << "LookAt: \"up\" vector (" << up.x << ", " << up.y << ", " << up.z
                << ") and viewing direction (" << dir.x << ", " << dir.y << ", " << dir.z
                << " passed to LookAt are pointing in the same direction." << std::endl;
        TanVector3f right = Normalize(Cross(Normalize(up), dir));
        TanVector3f newUp = Cross(dir, right);
        worldFromCamera[0][0] = right.x;
        worldFromCamera[1][0] = right.y;
        worldFromCamera[2][0] = right.z;
        worldFromCamera[3][0] = 0.f;
        worldFromCamera[0][1] = newUp.x;
        worldFromCamera[1][1] = newUp.y;
        worldFromCamera[2][1] = newUp.z;
        worldFromCamera[3][1] = 0.f;
        worldFromCamera[0][2] = dir.x;
        worldFromCamera[1][2] = dir.y;
        worldFromCamera[2][2] = dir.z;
        worldFromCamera[3][2] = 0.f;

        SquareMatrix<4> cameraFromWorld = InvertOrExit(worldFromCamera);
        return Transform(cameraFromWorld, worldFromCamera);
    }

    Transform Orthographic(float zNear, float zFar) {
        return Scale(1.f, 1.f, 1 / (zFar - zNear)) * Translate(TanVector3f(0.f, 0.f, -zNear));
    }

    Transform Perspective(float fov, float n, float f) {
        // Perform projective divide for perspective projection
        SquareMatrix<4> persp(1.f, 0.f, 0.f, 0.f,
            0.f, 1.f, 0.f, 0.f,
            0.f, 0.f, f / (f - n), -f * n / (f - n),
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

    void Transform::Decompose(TanVector3f* T, SquareMatrix<4>* R, SquareMatrix<4>* S) const {
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

} // namespace lightfold
