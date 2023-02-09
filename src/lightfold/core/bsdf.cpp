#include <core/bsdf.h>

namespace lightfold {

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