#pragma once

#include <core/bsdf.h>
#include <math/bsdfmath.h>

#include <iostream>

namespace lightfold {

    class Diffuse : public BSDF {
    public:
        Diffuse(const SurfaceInteraction& si, Spectrum color, float roughness) :
            color(color), roughness(roughness), BSDF(si, false) {}
        Spectrum distF(const Tangent3f& wo, const Tangent3f& wi) const;
        Spectrum Sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample, float* pdf, bool isSingular = false) const;
        //Spectrum reflectance(const Tangent3f& wo, int nSamples, const Point2f* samples) const;
        //Spectrum reflectance(int nSamples, const Point2f* samples1, const Point2f* samples2) const;
        float Pdf(const Tangent3f& wo, const Tangent3f& wi) const;

    private:
        // Diffuse Private Data
        Spectrum color;
        float roughness;
    };

    Spectrum Diffuse::distF(const Tangent3f& wo, const Tangent3f& wi) const {
        if (Dot(ng, wi) <= 0) {
            return { 0 };
        }
        const float FV = schlick_fresnel(Dot(ns, wo));
        const float FL = schlick_fresnel(Dot(ns, wi));
        float f = 0.0f;

        // Lambertian component.
        f += (1.0f - 0.5f * FV) * (1.0f - 0.5f * FL);

        // Retro-reflection component.
        const float LH2 = Dot(wi, wo) + 1;
        const float RR = roughness * LH2;
        f += RR * (FL + FV + FL * FV * (RR - 1.0f));

        float value = InvPi * Dot(ns, wi) * f;
        return color * Spectrum(value);
    }

    Spectrum Diffuse::Sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample, float* pdf, bool isSingular) const {
        CosineSampleHemisphere(sample);
        if (wo.z > 0) {
            *pdf = Pdf(wo, *wi);
            return distF(wo, *wi);
        }
        else {
            *pdf = 0.0f;
            return { 0 };
        }
    }

    float Diffuse::Pdf(const Tangent3f& wo, const Tangent3f& wi) const {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
    }

} // namespace lightfold