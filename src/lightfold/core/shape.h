#pragma once

#include <core/interaction.h>
#include <math/transform.h>

#include <vector>
#include <memory>

namespace lightfold {

    class Shape {
    public:
        // Shape Public Methods
        Shape(const Transform* ObjectToWorld, const Transform* WorldToObject, bool reverseOrientation) :
            ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject),
            reverseOrientation(reverseOrientation), transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}
        virtual Bounds3f ObjectBound() const = 0;
        virtual Bounds3f WorldBound() const;
        virtual bool Intersect(const Ray& ray, float* tHit, SurfaceInteraction* isect) const = 0;
        virtual bool IntersectP(const Ray& ray) const {
            float tHit = ray.tMax;
            SurfaceInteraction isect;
            return Intersect(ray, &tHit, &isect);
        }

        // Shape Public Data
        const Transform* ObjectToWorld, * WorldToObject;
        const bool reverseOrientation;
        const bool transformSwapsHandedness;
    };

    // Shape Methods Definitions
    inline Bounds3f Shape::WorldBound() const {
        return (*ObjectToWorld)(ObjectBound());
    }

} // namespace lightfold