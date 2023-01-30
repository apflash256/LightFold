#pragma once

#include <core/interaction.h>
#include <core/spectrum.h>

#include <math/sample.h>

namespace lightfold {

    class BSDF {
    public:
        // BSDF Interface
        virtual ~BSDF() {}
        BSDF() {}

        virtual Spectrum distF(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        // Gives the value of the distribution function for the given pair of directions
        virtual Spectrum Sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample,
            float* pdf) const = 0;
        // Samples an incident direction from the given outgoing direction and returns the value of the BxDF.
        virtual Spectrum reflectance(const Tangent3f& wo, int nSamples, const Point2f* samples) const;
        // Compute hemispherical-directional reflectance
        virtual Spectrum reflectance(int nSamples, const Point2f* samples1, const Point2f* samples2) const;
        // Compute hemispherical-hemispherical reflectance
        virtual float Pdf(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        // Gives the value of the PDF for the given pair of directions
    };

    Spectrum BSDF::reflectance(const Tangent3f& w, int nSamples, const Point2f* u) const {
        Spectrum r(0.f);
        for (int i = 0; i < nSamples; ++i) {
            // Estimate one term of $\rho_\roman{hd}$
            Tangent3f wi;
            float pdf = 0;
            Spectrum f = Sample_F(w, &wi, u[i], &pdf);
            if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
        }
        return r / nSamples;
    }

    Spectrum BSDF::reflectance(int nSamples, const Point2f* u1, const Point2f* u2) const {
        Spectrum r(0.f);
        for (int i = 0; i < nSamples; ++i) {
            // Estimate one term of $\rho_\roman{hh}$
            Tangent3f wo, wi;
            wo = UniformSampleHemisphere(u1[i]);
            float pdfo = UniformHemispherePdf(), pdfi = 0;
            Spectrum f = Sample_F(wo, &wi, u2[i], &pdfi);
            if (pdfi > 0)
                r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
        }
        return r / (Pi * nSamples);
    }

} // namespace lightfold