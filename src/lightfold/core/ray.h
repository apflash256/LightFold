#pragma once

#include <math/vecmath.h>

namespace lightfold {

    // Ray Definition
    class Ray {
    public:
        // Ray Public Methods
        Point3f operator()(float t) const { return o + d * t; }

        Ray() = default;
        Ray(Point3f o, TanVector3f d) : o(o), d(d) {}

        // Ray Public Members
        Point3f o;
        TanVector3f d;
    };

    // RayDifferential Definition
    class RayDifferential : public Ray {
    public:
        // RayDifferential Public Methods
        RayDifferential() = default;
        RayDifferential(Point3f o, TanVector3f d) : Ray(o, d) {}

        explicit RayDifferential(const Ray& ray) : Ray(ray) {}

        void ScaleDifferentials(float s) {
            rxOrigin = o + (rxOrigin - o) * s;
            ryOrigin = o + (ryOrigin - o) * s;
            rxDirection = d + (rxDirection - d) * s;
            ryDirection = d + (ryDirection - d) * s;
        }

        // RayDifferential Public Members
        bool hasDifferentials = false;
        Point3f rxOrigin, ryOrigin;
        TanVector3f rxDirection, ryDirection;
    };

    // Ray Inline Functions
    inline Point3f OffsetRayOrigin(Point3fi pi, CotVector3f n, TanVector3f w) {
        // Find vector _offset_ to corner of error bounds and compute initial _po_
        float d = Dot(Abs(n), pi.Error());
        TanVector3f offset = d * TanVector3f(n);
        if (Dot(w, n) < 0)
            offset = -offset;
        Point3f po = Point3f(pi) + offset;
        return po;
    }

    inline Ray SpawnRay(Point3fi pi, CotVector3f n, TanVector3f d) {
        return Ray(OffsetRayOrigin(pi, n, d), d);
    }

    inline Ray SpawnRayTo(Point3fi pFrom, CotVector3f n, Point3f pTo) {
        TanVector3f d = pTo - Point3f(pFrom);
        return SpawnRay(pFrom, n, d);
    }

    inline Ray SpawnRayTo(Point3fi pFrom, CotVector3f nFrom, Point3fi pTo, CotVector3f nTo) {
        Point3f pf = OffsetRayOrigin(pFrom, nFrom, Point3f(pTo) - Point3f(pFrom));
        Point3f pt = OffsetRayOrigin(pTo, nTo, pf - Point3f(pTo));
        return Ray(pf, pt - pf);
    }

} // namespace lightfold