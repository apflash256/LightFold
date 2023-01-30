#include <core/interaction.h>
#include <core/light.h>

namespace lightfold {

    SurfaceInteraction::SurfaceInteraction(
        const Point3f& p, const Tangent3f& pError, const Point2f& uv,
        const Tangent3f& wo, const Tangent3f& dpdu, const Tangent3f& dpdv,
        const Normal3f& dndu, const Normal3f& dndv, float time, const Shape* shape)
        : Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), pError, wo, time,
            nullptr),
        uv(uv),
        dpdu(dpdu),
        dpdv(dpdv),
        dndu(dndu),
        dndv(dndv),
        shape(shape) {
        // Initialize shading geometry from true geometry
        shading.n = n;
        shading.dpdu = dpdu;
        shading.dpdv = dpdv;
        shading.dndu = dndu;
        shading.dndv = dndv;

        // Adjust normal based on orientation and handedness
        if (shape &&
            (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
            n *= -1;
            shading.n *= -1;
        }
    }

    void SurfaceInteraction::SetShadingGeometry(const Tangent3f& dpdus,
        const Tangent3f& dpdvs,
        const Normal3f& dndus,
        const Normal3f& dndvs,
        bool orientationIsAuthoritative) {
        // Compute _shading.n_ for _SurfaceInteraction_
        shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
        if (orientationIsAuthoritative)
            n = FaceForward(n, shading.n);
        else
            shading.n = FaceForward(shading.n, n);

        // Initialize _shading_ partial derivative values
        shading.dpdu = dpdus;
        shading.dpdv = dpdvs;
        shading.dndu = dndus;
        shading.dndv = dndvs;
    }

    void SurfaceInteraction::ComputeScatteringFunctions(const RayDifferential& ray,
        bool allowMultipleLobes, TransportMode mode) {
        ComputeDifferentials(ray);
        primitive->ComputeScatteringFunctions(this, mode,
            allowMultipleLobes);
    }

    void SurfaceInteraction::ComputeDifferentials(
        const RayDifferential& ray) const {
        if (ray.hasDifferentials) {
            // Estimate screen space change in $\pt{}$ and $(u,v)$

            // Compute auxiliary intersection points with plane
            float tx =
                -Dot(n, (ray.rxOrigin - p)) / Dot(n, ray.rxDirection);
            if (std::isinf(tx) || std::isnan(tx)) goto fail;
            Point3f px = ray.rxOrigin + tx * ray.rxDirection;
            float ty =
                -Dot(n, (ray.ryOrigin - p)) / Dot(n, ray.ryDirection);
            if (std::isinf(ty) || std::isnan(ty)) goto fail;
            Point3f py = ray.ryOrigin + ty * ray.ryDirection;
            dpdx = px - p;
            dpdy = py - p;

            // Compute $(u,v)$ offsets at auxiliary points

            // Choose two dimensions to use for ray offset computation
            int dim[2];
            if (std::abs(n.x) > std::abs(n.y) && std::abs(n.x) > std::abs(n.z)) {
                dim[0] = 1;
                dim[1] = 2;
            }
            else if (std::abs(n.y) > std::abs(n.z)) {
                dim[0] = 0;
                dim[1] = 2;
            }
            else {
                dim[0] = 0;
                dim[1] = 1;
            }

            // Initialize _A_, _Bx_, and _By_ matrices for offset computation
            float A[2][2] = { {dpdu[dim[0]], dpdv[dim[0]]},
                             {dpdu[dim[1]], dpdv[dim[1]]} };
            float Bx[2] = { px[dim[0]] - p[dim[0]], px[dim[1]] - p[dim[1]] };
            float By[2] = { py[dim[0]] - p[dim[0]], py[dim[1]] - p[dim[1]] };
            if (!SolveLinearSystem2x2(A, Bx, &dudx, &dvdx)) dudx = dvdx = 0;
            if (!SolveLinearSystem2x2(A, By, &dudy, &dvdy)) dudy = dvdy = 0;
        }
        else {
        fail:
            dudx = dvdx = 0;
            dudy = dvdy = 0;
            dpdx = dpdy = Tangent3f(0, 0, 0);
        }
    }

    Spectrum SurfaceInteraction::Le(const Tangent3f& w) const {
        const AreaLight* area = primitive->GetAreaLight();
        return area ? area->L(*this, w) : Spectrum(0.f);
    }

} // namespace lightfold