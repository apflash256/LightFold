#pragma once

#include <math/vecmath.h>

namespace lightfold {

    class Medium;

    class Ray {
    public:
        // Ray Public Methods
        Ray() : tMax(Infinity), time(0.f), medium(nullptr) { }
        Ray(const Point3f& o, const Tangent3f& d, float tMax = Infinity,
            float time = 0.f, const Medium* medium = nullptr)
            : o(o), d(d), tMax(tMax), time(time), medium(medium) { }
        Point3f operator()(float t) const { return o + d * t; }

        // Ray Public Data
        Point3f o;
        Tangent3f d;
        mutable float tMax;
        float time;
        const Medium* medium;
    };

    class RayDifferential : public Ray {
    public:
        // RayDifferential Public Methods
        RayDifferential() { hasDifferentials = false; }
        RayDifferential(const Point3f& o, const Tangent3f& d,
            float tMax = Infinity, float time = 0.f,
            const Medium* medium = nullptr)
            : Ray(o, d, tMax, time, medium) {
            hasDifferentials = false;
        }
        RayDifferential(const Ray& ray) : Ray(ray) {
            hasDifferentials = false;
        }
        void ScaleDifferentials(float s) {
            rxOrigin = o + (rxOrigin - o) * s;
            ryOrigin = o + (ryOrigin - o) * s;
            rxDirection = d + (rxDirection - d) * s;
            ryDirection = d + (ryDirection - d) * s;
        }

        // RayDifferential Public Data
        bool hasDifferentials;
        Point3f rxOrigin, ryOrigin;
        Tangent3f rxDirection, ryDirection;

    };

} // namespace lightfold