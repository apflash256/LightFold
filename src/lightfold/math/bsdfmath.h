#pragma once

#include <math/vecmath.h>
#include <core/spectrum.h>

namespace lightfold {

    inline Tangent3f Reflect(const Tangent3f& wo, const Tangent3f& n) {
        return -wo + 2 * Dot(wo, n) * n;
    }

    inline bool Refract(const Tangent3f& wi, const Normal3f& n, float eta,
        Tangent3f* wt) {
        // Compute $\cos \theta_\roman{t}$ using Snell's law
        float cosThetaI = Dot(n, wi);
        float sin2ThetaI = std::max(float(0), float(1 - cosThetaI * cosThetaI));
        float sin2ThetaT = eta * eta * sin2ThetaI;

        // Handle total internal reflection for transmission
        if (sin2ThetaT >= 1) return false;
        float cosThetaT = std::sqrt(1 - sin2ThetaT);
        *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Tangent3f(n);
        return true;
    }

    float fresnel_dielectric(float eta, const Normal3f N, const Tangent3f I,
        Tangent3f* R, Tangent3f* T, bool* is_inside)
    {
        float cos = Dot(N, I), neta;
        Normal3f Nn;

        // check which side of the surface we are on
        if (cos > 0) {
            // we are on the outside of the surface, going in
            neta = 1 / eta;
            Nn = N;
            *is_inside = false;
        }
        else {
            // we are inside the surface
            cos = -cos;
            neta = eta;
            Nn = -N;
            *is_inside = true;
        }

        // compute reflection
        *R = (2 * cos) * (Tangent3f)Nn - I;

        float arg = 1 - (neta * neta * (1 - (cos * cos)));
        if (arg < 0) {
            *T = Tangent3f(0.0f, 0.0f, 0.0f);
            return 1;  // total internal reflection
        }
        else {
            float dnp = std::max(sqrtf(arg), 1e-7f);
            float nK = (neta * cos) - dnp;
            *T = -(neta * I) + (nK * (Tangent3f)Nn);
            // compute Fresnel terms
            float cosTheta1 = cos;  // N.R
            float cosTheta2 = -Dot(Nn, *T);
            float pPara = (cosTheta1 - eta * cosTheta2) / (cosTheta1 + eta * cosTheta2);
            float pPerp = (eta * cosTheta1 - cosTheta2) / (eta * cosTheta1 + cosTheta2);
            return 0.5f * (pPara * pPara + pPerp * pPerp);
        }
    }

    float fresnel_dielectric_cos(float cosi, float eta)
    {
        // compute fresnel reflectance without explicitly computing
        // the refracted direction
        float c = fabsf(cosi);
        float g = eta * eta - 1 + c * c;
        if (g > 0) {
            g = sqrtf(g);
            float A = (g - c) / (g + c);
            float B = (c * (g + c) - 1) / (c * (g - c) + 1);
            return 0.5f * A * A * (1 + B * B);
        }
        return 1.0f;  // TIR(no refracted component)
    }

    Spectrum fresnel_conductor(float cosi, const Spectrum eta, const Spectrum k)
    {
        Spectrum tmp_f = eta * eta + k * k;
        Spectrum tmp = tmp_f * cosi;
        Spectrum Rparl2 = (tmp - (2.0f * eta * cosi) + 1.0f) / (tmp + (2.0f * eta * cosi) + 1.0f);
        Spectrum Rperp2 = (tmp_f - (2.0f * eta * cosi) + cosi) / (tmp_f + (2.0f * eta * cosi) + cosi);
        return (Rparl2 + Rperp2) * 0.5f;
    }

    float schlick_fresnel(float u)
    {
        float m = Clamp(1.0f - u, 0.0f, 1.0f);
        float m2 = m * m;
        return m2 * m2 * m;  // pow(m, 5)
    }

} // namespace lightfold