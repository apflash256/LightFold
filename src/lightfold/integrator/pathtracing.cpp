#include <integrator/pathtracing.h>

#include <core/bsdf.h>
#include <core/stats.h>

#include <iostream>

namespace lightfold {

    Spectrum PathIntegrator::Li(const RayDifferential& r, const Scene& scene,
        Sampler& sampler, int depth) const {
        Spectrum L(0.f), beta(1.f);
        RayDifferential ray(r);
        bool specularBounce = false;
        for (int bounces = 0; ; ++bounces) {
            SurfaceInteraction isect;
            bool foundIntersection = scene.Intersect(ray, &isect);
            if (bounces == 0 || specularBounce) {
                if (foundIntersection)
                    L += beta * isect.Le(-ray.d);
                else
                    for (const auto& light : scene.lights)
                        L += beta * light->Le(ray);
            }
            if (!foundIntersection || bounces >= maxDepth)
                break;
            isect.ComputeScatteringFunctions(ray, true);
            if (!isect.bsdf) {
                ray = isect.SpawnRay(ray.d);
                bounces--;
                continue;
            }
            L += beta * UniformSampleOneLight(isect, scene, sampler);
            Tangent3f wo = -ray.d, wi;
            float pdf;
            Spectrum f = isect.bsdf->Sample_F(wo, &wi, sampler.Get2D(), &pdf);
            if (f.IsBlack() || pdf == 0.f) {
                free(isect.bsdf);
                break;
            }
            beta *= f * AbsDot(wi, isect.shading.n) / pdf;
            ray = isect.SpawnRay(wi);
            // Account for subsurface scattering, if applicable
            // Something?
            if (bounces > 3) {
                float q = std::max(.05f, 1 - beta.Value());
                if (sampler.Get1D() < q) {
                    free(isect.bsdf);
                    break;
                }
                beta /= 1 - q;
            }
            free(isect.bsdf);
        }
        return L;
    }

} // namespace lightfold