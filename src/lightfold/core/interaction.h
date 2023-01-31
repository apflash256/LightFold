#pragma once

#include <core/material.h>
#include <core/medium.h>

namespace lightfold {

    class Shape;
    class Primitive;

    struct Interaction {
        // Interaction Public Methods
        Interaction() : time(0) { }
        Interaction(const Point3f& p, const Normal3f& n, const Tangent3f& pError,
            const Tangent3f& wo, float time,
            const MediumInterface& mediumInterface)
            : p(p), time(time), pError(pError), wo(wo), n(n),
            mediumInterface(mediumInterface) { }
        bool IsSurfaceInteraction() const {
            return n != Normal3f();
        }
        Ray SpawnRay(const Tangent3f& d) const {
            Point3f o = OffsetRayOrigin(p, pError, n, d);
            return Ray(o, d, Infinity, time, GetMedium(d));
        }
        Ray SpawnRayTo(const Point3f& p2) const {
            Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
            Tangent3f d = p2 - origin;
            return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
        }
        Ray SpawnRayTo(const Interaction& it) const {
            Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
            Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
            Tangent3f d = target - origin;
            return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
        }
        Interaction(const Point3f& p, const Tangent3f& wo, float time,
            const MediumInterface& mediumInterface)
            : p(p), time(time), wo(wo), mediumInterface(mediumInterface) { }
        Interaction(const Point3f& p, float time,
            const MediumInterface& mediumInterface)
            : p(p), time(time), mediumInterface(mediumInterface) { }
        bool IsMediumInteraction() const { return !IsSurfaceInteraction(); }
        const Medium* GetMedium(const Tangent3f& w) const {
            return Dot(w, n) > 0 ? mediumInterface.outside :
                mediumInterface.inside;
        }
        const Medium* GetMedium() const {
            return mediumInterface.inside;
        }

        // Interaction Public Data
        Point3f p;
        float time;
        Tangent3f pError;
        Tangent3f wo;
        Normal3f n;
        MediumInterface mediumInterface;
    };

    class Material;
    class BSDF;
    class BSSRDF;

    class SurfaceInteraction : public Interaction {
    public:
        // SurfaceInteraction Public Methods
        SurfaceInteraction() { }
        SurfaceInteraction(const Point3f& p, const Tangent3f& pError, const Point2f& uv,
            const Tangent3f& wo, const Tangent3f& dpdu, const Tangent3f& dpdv,
            const Normal3f& dndu, const Normal3f& dndv,
            float time, const Shape* sh);
        void SetShadingGeometry(const Tangent3f& dpdu, const Tangent3f& dpdv,
            const Normal3f& dndu, const Normal3f& dndv, bool orientationIsAuthoritative);
        void ComputeScatteringFunctions(const RayDifferential& ray,
            bool allowMultipleLobes = false, TransportMode mode = TransportMode::Radiance);
        void ComputeDifferentials(const RayDifferential& r) const;
        Spectrum Le(const Tangent3f& w) const;

        // SurfaceInteraction Public Data
        Point2f uv;
        Tangent3f dpdu, dpdv;
        Normal3f dndu, dndv;
        const Shape* shape = nullptr;
        struct {
            Normal3f n;
            Tangent3f dpdu, dpdv;
            Normal3f dndu, dndv;
        } shading; // Shading Geometry
        const Primitive* primitive = nullptr;
        BSDF* bsdf = nullptr;
        BSSRDF* bssrdf = nullptr;
        mutable Tangent3f dpdx, dpdy;
        mutable float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;
    };

    class MediumInteraction : public Interaction {
    public:
        // MediumInteraction Public Methods
        MediumInteraction(const Point3f& p, const Tangent3f& wo, float time,
            const Medium* medium, const PhaseFunction* phase)
            : Interaction(p, wo, time, medium), phase(phase) { }
        bool IsValid() const { return phase != nullptr; }

        // MediumInteraction Public Data
        const PhaseFunction* phase;
    };

} // namespace lightfold