#include <core/scene.h>


#include <core/interaction.h>
#include <core/medium.h>
#include <core/sampler.h>

namespace lightfold {

    bool Scene::Intersect(const Ray& ray, SurfaceInteraction* isect) const {
        return aggregate->Intersect(ray, isect);
    }

    bool Scene::IntersectP(const Ray& ray) const {
        return aggregate->IntersectP(ray);
    }

    bool Scene::IntersectTr(Ray ray, Sampler& sampler, SurfaceInteraction* isect,
        Spectrum* Tr) const {
        *Tr = Spectrum(1.f);
        while (true) {
            bool hitSurface = Intersect(ray, isect);
            // Accumulate beam transmittance for ray segment
            if (ray.medium) *Tr *= ray.medium->Transmittance(ray, sampler);

            // Initialize next ray segment or terminate transmittance computation
            if (!hitSurface) return false;
            if (isect->primitive->GetMaterial() != nullptr) return true;
            ray = isect->SpawnRay(ray.d);
        }
    }

} // namespace lightfold