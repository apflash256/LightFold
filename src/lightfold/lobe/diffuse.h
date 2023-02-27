#pragma once

#include <core/lobe.h>
#include <math/bsdfmath.h>

#include <iostream>

namespace lightfold {

    class Diffuse : public Lobe {
    public:
        Diffuse(Spectrum color, float roughness) : Lobe(LobeType(LOBE_REFLECTION | LOBE_DIFFUSE)),
            color(color), roughness(roughness) {}
        Spectrum dist_F(const Tangent3f& wo, const Tangent3f& wi) const;
        Spectrum sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample,
            float* pdf, LobeType* sampledType = nullptr) const;
        float PDF(const Tangent3f& wo, const Tangent3f& wi) const;

    private:
        // Diffuse Private Data
        Spectrum color;
        float roughness;
    };

    Spectrum Diffuse::dist_F(const Tangent3f& wo, const Tangent3f& wi) const {
        const float FV = schlick_fresnel(CosTheta(wo));
        const float FL = schlick_fresnel(CosTheta(wi));
        float f = 0.0f;

        // Lambertian component.
        f += (1.0f - 0.5f * FV) * (1.0f - 0.5f * FL);

        // Retro-reflection component.
        const float LH2 = Dot(wi, wo) + 1;
        const float RR = roughness * LH2;
        f += RR * (FL + FV + FL * FV * (RR - 1.0f));

        float value = InvPi * (CosTheta(wi)) * f;
        return color * Spectrum(value);
    }

    Spectrum Diffuse::sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample,
        float* pdf, LobeType* sampledType) const {
        *wi = CosineSampleHemisphere(sample);
        if (wo.z < 0) wi->z *= -1;
        *pdf = PDF(wo, *wi);
        return dist_F(wo, *wi);
    }

    float Diffuse::PDF(const Tangent3f& wo, const Tangent3f& wi) const {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
    }

} // namespace lightfold