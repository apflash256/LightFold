#pragma once

#include <core/material.h>
#include <core/spectrum.h>

namespace lightfold {

    class MatteMaterial : public Material {
    public:
        // MatteMaterial Public Methods
        MatteMaterial(const std::shared_ptr<Texture<Spectrum>>& color,
            const std::shared_ptr<Texture<float>>& roughness)
            : color(color), roughness(roughness) {}
        void ComputeScatteringFunctions(SurfaceInteraction* si, TransportMode mode,
            bool allowMultipleLobes) const;

    private:
        // MatteMaterial Private Data
        std::shared_ptr<Texture<Spectrum>> color;
        std::shared_ptr<Texture<float>> roughness;
    };

} // namespace lightfold