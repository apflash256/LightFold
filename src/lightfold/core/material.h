#pragma once

#include <core/texture.h>

#include <memory>

namespace lightfold {

	enum class TransportMode { Radiance, Importance };

    class Material {
    public:
        // Material Interface
        virtual void ComputeScatteringFunctions(SurfaceInteraction* si, TransportMode mode,
            bool allowMultipleLobes) const = 0;
        virtual ~Material() {}
        static void Bump(const std::shared_ptr<Texture<float>>& d, SurfaceInteraction* si);
    };

} // namespace lightfold