#include <objects/shape.h>

namespace lightfold {

    std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
        const Transform* ObjectToWorld, const Transform* WorldToObject,
        bool reverseOrientation, int nTriangles,
        const int* vertexIndices, int nVertices, const Point3f* p,
        const TanVector3f* s, const CotVector3f* n, const Point2f* uv) {
        std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
            *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv);
        std::vector<std::shared_ptr<Shape>> tris;
        for (int i = 0; i < nTriangles; ++i)
            tris.push_back(std::make_shared<Triangle>(ObjectToWorld,
                WorldToObject, reverseOrientation, mesh, i));
        return tris;
    }

    Bounds3f Triangle::ObjectBound() const {
            const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
            (*WorldToObject)(p2));
    }

    Bounds3f Triangle::WorldBound() const {
            const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        return Union(Bounds3f(p0, p1), p2);
    }

    bool Triangle::Intersect(const Ray &ray, float *tHit) const {
        const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];
        Point3f p0t = p0 - (TanVector3f)(GeoVector3f)ray.o;
        Point3f p1t = p1 - (TanVector3f)(GeoVector3f)ray.o;
        Point3f p2t = p2 - (TanVector3f)(GeoVector3f)ray.o;
        int kz = MaxComponentIndex(Abs(ray.d));
        int kx = kz + 1; if (kx == 3) kx = 0;
        int ky = kx + 1; if (ky == 3) ky = 0;
        TanVector3f d = Permute(ray.d, { kx, ky, kz });
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
        //float b0 = e0 * invDet;
        //float b1 = e1 * invDet;
        //float b2 = e2 * invDet;
        float t = tScaled * invDet;
        float maxZt = MaxComponentValue(Abs(TanVector3f(p0t.z, p1t.z, p2t.z)));
        float deltaZ = gamma(3) * maxZt;
        float maxXt = MaxComponentValue(Abs(TanVector3f(p0t.x, p1t.x, p2t.x)));
        float maxYt = MaxComponentValue(Abs(TanVector3f(p0t.y, p1t.y, p2t.y)));
        float deltaX = gamma(5) * (maxXt + maxZt);
        float deltaY = gamma(5) * (maxYt + maxZt);
        float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
        float maxE = MaxComponentValue(Abs(TanVector3f(e0, e1, e2)));
        float deltaT = 3 * (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * std::abs(invDet);
        if (t <= deltaT)
            return false;
        *tHit = t;
        return true;
}

} // namespace lightfold