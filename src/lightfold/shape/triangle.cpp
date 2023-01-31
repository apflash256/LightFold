#include <shape/triangle.h>

namespace lightfold {


    std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
        const Transform* ObjectToWorld, const Transform* WorldToObject,
        bool reverseOrientation, int nTriangles, const int* vertexIndices,
        int nVertices, const Point3f* p, const Tangent3f* s, const Normal3f* n,
        const Point2f* uv) {
        std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
            *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv);
        std::vector<std::shared_ptr<Shape>> tris;
        tris.reserve(nTriangles);
        for (int i = 0; i < nTriangles; ++i)
            tris.push_back(std::make_shared<Triangle>(ObjectToWorld,
                WorldToObject, reverseOrientation, mesh, i));
        return tris;
    }

    Bounds3f Triangle::ObjectBound() const {
        const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)), (*WorldToObject)(p2));
    }

    Bounds3f Triangle::WorldBound() const {
        const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        return Union(Bounds3f(p0, p1), p2);
    }

    bool Triangle::Intersect(const Ray& ray, float* tHit, SurfaceInteraction* isect) const {
        const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        Point3f p0t = p0 - (Tangent3f)(Vector3f)ray.o;
        Point3f p1t = p1 - (Tangent3f)(Vector3f)ray.o;
        Point3f p2t = p2 - (Tangent3f)(Vector3f)ray.o;
        int kz = MaxComponentIndex(Abs(ray.d));
        int kx = kz + 1; if (kx == 3) kx = 0;
        int ky = kx + 1; if (ky == 3) ky = 0;
        Tangent3f d = Permute(ray.d, { kx, ky, kz });
        p0t = Permute(p0t, { kx, ky, kz });
        p1t = Permute(p1t, { kx, ky, kz });
        p2t = Permute(p2t, { kx, ky, kz });
        float Sx = -d.x / d.z;
        float Sy = -d.y / d.z;
        float Sz = 1.f / d.z;
        p0t.x += Sx * p0t.z;
        p0t.y += Sy * p0t.z;
        p1t.x += Sx * p1t.z;
        p1t.y += Sy * p1t.z;
        p2t.x += Sx * p2t.z;
        p2t.y += Sy * p2t.z;
        float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        float e2 = p0t.x * p1t.y - p0t.y * p1t.x;
        if (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f) {
            double p2txp1ty = (double)p2t.x * (double)p1t.y;
            double p2typ1tx = (double)p2t.y * (double)p1t.x;
            e0 = (float)(p2typ1tx - p2txp1ty);
            double p0txp2ty = (double)p0t.x * (double)p2t.y;
            double p0typ2tx = (double)p0t.y * (double)p2t.x;
            e1 = (float)(p0typ2tx - p0txp2ty);
            double p1txp0ty = (double)p1t.x * (double)p0t.y;
            double p1typ0tx = (double)p1t.y * (double)p0t.x;
            e2 = (float)(p1typ0tx - p1txp0ty);
        }
        if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
            return false;
        float det = e0 + e1 + e2;
        if (det == 0)
            return false;
        p0t.z *= Sz;
        p1t.z *= Sz;
        p2t.z *= Sz;
        float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
            return false;
        else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
            return false;
        float invDet = 1 / det;
        float b0 = e0 * invDet;
        float b1 = e1 * invDet;
        float b2 = e2 * invDet;
        float t = tScaled * invDet;
        float maxZt = MaxComponentValue(Abs(Tangent3f(p0t.z, p1t.z, p2t.z)));
        float deltaZ = gamma(3) * maxZt;
        float maxXt = MaxComponentValue(Abs(Tangent3f(p0t.x, p1t.x, p2t.x)));
        float maxYt = MaxComponentValue(Abs(Tangent3f(p0t.y, p1t.y, p2t.y)));
        float deltaX = gamma(5) * (maxXt + maxZt);
        float deltaY = gamma(5) * (maxYt + maxZt);
        float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
        float maxE = MaxComponentValue(Abs(Tangent3f(e0, e1, e2)));
        float deltaT = 3 * (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * std::abs(invDet);
        if (t <= deltaT)
            return false;
        Tangent3f dpdu, dpdv;
        Point2f uv[3];
        GetUVs(uv);
        Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
        Tangent3f dp02 = p0 - p2, dp12 = p1 - p2;
        float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        if (determinant == 0) {
            CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);
        }
        else {
            float invdet = 1 / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        float xAbsSum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
        float yAbsSum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
        float zAbsSum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
        Tangent3f pError = gamma(7) * Tangent3f(xAbsSum, yAbsSum, zAbsSum);
        Point3f pHit(b0 * p0.x + b1 * p1.x + b2 * p2.x,
            b0 * p0.y + b1 * p1.y + b2 * p2.y,
            b0 * p0.z + b1 * p1.z + b2 * p2.z);
        Point2f uvHit(b0 * uv[0].x + b1 * uv[1].x + b2 * uv[2].x,
            b0 * uv[0].y + b1 * uv[1].y + b2 * uv[2].y);
        *isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
            Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);
        isect->n = isect->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
        if (mesh->n || mesh->s) {
            Normal3f ns;
            if (mesh->n) ns = Normalize(b0 * mesh->n[v[0]] +
                b1 * mesh->n[v[1]] +
                b2 * mesh->n[v[2]]);
            else
                ns = isect->n;
            Tangent3f ss;
            if (mesh->s) ss = Normalize(b0 * mesh->s[v[0]] +
                b1 * mesh->s[v[1]] +
                b2 * mesh->s[v[2]]);
            else
                ss = Normalize(isect->dpdu);
            Tangent3f ts = Cross(ns, ss);
            if (LengthSquared(ts) > 0.f) {
                ts = Normalize(ts);
                ss = Cross(ts, ns);
            }
            else
                CoordinateSystem((Tangent3f)ns, &ss, &ts);
            Normal3f dndu, dndv;
            if (mesh->n) {
                Vector2f duv02 = uv[0] - uv[2];
                Vector2f duv12 = uv[1] - uv[2];
                Normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
                Normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
                float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
                if (determinant == 0)
                    dndu = dndv = Normal3f(0, 0, 0);
                else {
                    float invDet = 1 / determinant;
                    dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
                    dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
                }
            }
            else
                dndu = dndv = Normal3f(0, 0, 0);
            isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
        }
        if (mesh->n)
            isect->n = FaceForward(isect->n, isect->shading.n);
        else if (reverseOrientation ^ transformSwapsHandedness)
            isect->n = isect->shading.n = -isect->n;
        *tHit = t;
        return true;
    }

} // namespace lightfold