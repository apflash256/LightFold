#pragma once

#include <core/integrator.h>

namespace lightfold {

    enum class LightStrategy { UniformSampleAll, UniformSampleOne };

    class DirectLightingIntegrator : public SamplerIntegrator {
    public:
        // DirectLightingIntegrator Public Methods
        DirectLightingIntegrator(LightStrategy strategy, int maxDepth,
            std::shared_ptr<const Camera> camera,
            std::shared_ptr<Sampler> sampler,
            const Bounds2i& pixelBounds)
            : SamplerIntegrator(camera, sampler, pixelBounds),
            strategy(strategy),
            maxDepth(maxDepth) {}
        Spectrum Li(const RayDifferential& ray, const Scene& scene,
            Sampler& sampler, int depth) const;
        void Preprocess(const Scene& scene, Sampler& sampler);

    private:
        // DirectLightingIntegrator Private Data
        const LightStrategy strategy;
        const int maxDepth;
        std::vector<int> nLightSamples;
    };

} // namespace lightfold