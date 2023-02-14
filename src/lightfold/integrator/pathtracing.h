#pragma once

#include <core/integrator.h>

namespace lightfold {

    class PathIntegrator : public SamplerIntegrator {
    public:
        // PathIntegrator Public Methods
        Spectrum Li(const RayDifferential& ray, const Scene& scene, Sampler& sampler,
            int depth) const;
        PathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
            std::shared_ptr<Sampler> sampler, const Bounds2i& pixelBounds)
            : SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) { }

    private:
        // PathIntegrator Private Data
        const int maxDepth;
    };

} // namespace lightfold