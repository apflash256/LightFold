#pragma once

#include <core/material.h>
#include <math/ray.h>

namespace lightfold {

    class AreaLight;
    class SurfaceInteraction;

    class Primitive {
    public:
        // Primitive Interface
        virtual ~Primitive() { };
        virtual Bounds3f WorldBound() const = 0;
        virtual bool Intersect(const Ray& r, SurfaceInteraction*) const = 0;
        virtual bool IntersectP(const Ray& r) const = 0;
        virtual const AreaLight* GetAreaLight() const = 0;
        virtual const Material* GetMaterial() const = 0;
        virtual void ComputeScatteringFunctions(SurfaceInteraction* isect, TransportMode mode,
            bool allowMultipleLobes) const = 0;
    };

} // namespace lightfold