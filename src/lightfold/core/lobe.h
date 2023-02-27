#pragma once

#include <core/spectrum.h>
#include <math/sample.h>

namespace lightfold {

    enum LobeType {
        LOBE_REFLECTION = 1 << 0,
        LOBE_TRANSMISSION = 1 << 1,
        LOBE_DIFFUSE = 1 << 2,
        LOBE_GLOSSY = 1 << 3,
        LOBE_SPECULAR = 1 << 4,
        LOBE_ALL = LOBE_DIFFUSE | LOBE_GLOSSY | LOBE_SPECULAR | LOBE_REFLECTION | LOBE_TRANSMISSION,
    };

    class Lobe {
    public:
        // Lobe Interface
        virtual ~Lobe() {}
        Lobe(LobeType type) : type(type) {}
        bool matchesFlags(LobeType t) const { return (type & t) == type; }
        virtual Spectrum dist_F(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        virtual Spectrum sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample,
            float* pdf, LobeType* sampledType = nullptr) const = 0;
        virtual Spectrum reflectance(const Tangent3f& wo, int nSamples, const Point2f* samples) const;
        virtual Spectrum reflectance(int nSamples, const Point2f* samples1,
            const Point2f* samples2) const;
        virtual float PDF(const Tangent3f& wo, const Tangent3f& wi) const = 0;

        // BSDF Public Data
        const LobeType type;
    };


} // namespace lightfold