#include <core/primitive.h>
#include <core/stats.h>

namespace lightfold {

    const AreaLight* Aggregate::GetAreaLight() const {
        return nullptr;
    }

    const Material* Aggregate::GetMaterial() const {
        return nullptr;
    }

    void Aggregate::ComputeScatteringFunctions(SurfaceInteraction* isect,
        TransportMode mode, bool allowMultipleLobes) const {
    }

    bool TransformedPrimitive::Intersect(const Ray& r,
        SurfaceInteraction* isect) const {
        // Compute _ray_ after transformation by _PrimitiveToWorld_
        Transform InterpolatedPrimToWorld;
        PrimitiveToWorld.Interpolate(r.time, InterpolatedPrimToWorld);
        Ray ray = Inverse(InterpolatedPrimToWorld)(r);
        if (!primitive->Intersect(ray, isect)) return false;
        r.tMax = ray.tMax;
        // Transform instance's intersection data to world space
        if (!InterpolatedPrimToWorld.IsIdentity())
            *isect = InterpolatedPrimToWorld(*isect);
        return true;
    }

    bool TransformedPrimitive::IntersectP(const Ray& r) const {
        Transform InterpolatedPrimToWorld;
        PrimitiveToWorld.Interpolate(r.time, InterpolatedPrimToWorld);
        Transform InterpolatedWorldToPrim = Inverse(InterpolatedPrimToWorld);
        return primitive->IntersectP(InterpolatedWorldToPrim(r));
    }

    Bounds3f GeometricPrimitive::WorldBound() const { return shape->WorldBound(); }

    bool GeometricPrimitive::IntersectP(const Ray& r) const {
        return shape->IntersectP(r);
    }

    bool GeometricPrimitive::Intersect(const Ray& r,
        SurfaceInteraction* isect) const {
        float tHit;
        if (!shape->Intersect(r, &tHit, isect)) return false;
        r.tMax = tHit;
        isect->primitive = this;
        // Initialize _SurfaceInteraction::mediumInterface_ after _Shape_
        // intersection
        if (mediumInterface.IsMediumTransition())
            isect->mediumInterface = mediumInterface;
        else
            isect->mediumInterface = MediumInterface(r.medium);
        return true;
    }

    const AreaLight* GeometricPrimitive::GetAreaLight() const {
        return areaLight.get();
    }

    const Material* GeometricPrimitive::GetMaterial() const {
        return material.get();
    }

    void GeometricPrimitive::ComputeScatteringFunctions(
        SurfaceInteraction* isect, TransportMode mode, bool allowMultipleLobes) const {
        ProfilePhase p(Prof::ComputeScatteringFuncs);
        if (material)
            material->ComputeScatteringFunctions(isect, mode, allowMultipleLobes);
    }

} // namespace lightfold