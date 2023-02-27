#include <core/lobe.h>

#include <core/stats.h>

namespace lightfold {

    Spectrum Lobe::reflectance(const Tangent3f& w, int nSamples, const Point2f* u) const {
        Spectrum r(0.f);
        for (int i = 0; i < nSamples; ++i) {
            Tangent3f wi;
            float pdf = 0;
            Spectrum f = sample_F(w, &wi, u[i], &pdf);
            if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
        }
        return r / nSamples;
    }

    Spectrum Lobe::reflectance(int nSamples, const Point2f* u1, const Point2f* u2) const {
        Spectrum r(0.f);
        for (int i = 0; i < nSamples; ++i) {
            Tangent3f wo, wi;
            wo = UniformSampleHemisphere(u1[i]);
            float pdfo = UniformHemispherePdf(), pdfi = 0;
            Spectrum f = sample_F(wo, &wi, u2[i], &pdfi);
            if (pdfi > 0)
                r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
        }
        return r / (Pi * nSamples);
    }

} // namespace lightfold