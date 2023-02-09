#pragma once

#include <math/vecmath.h>

namespace lightfold {

    inline Point2f ConcentricSampleDisk(const Tuple2<float>& u) {
        Tuple2<float> uOffset(2 * u[0] - 1, 2 * u[1] - 1);
        if (uOffset.x == 0 && uOffset.y == 0)
            return Point2f(0, 0);
        float theta, r;
        if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
            r = uOffset.x;
            theta = PiOver4 * (uOffset.y / uOffset.x);
        }
        else {
            r = uOffset.y;
            theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
        }
        return Point2f(r*std::cos(theta), r*std::sin(theta));
    }

    inline Tangent3f UniformSampleHemisphere(const Tuple2<float>& u) {
        float z = u[0];
        float r = std::sqrt(std::max(0.f, 1.f - z * z));
        float phi = 2 * Pi * u[1];
        return Tangent3f(r * std::cos(phi), r * std::sin(phi), z);
    }

    inline float UniformHemispherePdf() { return Inv2Pi; }

    inline Tangent3f UniformSampleSphere(const Point2f& u) {
        float z = 1 - 2 * u[0];
        float r = std::sqrt(std::max(0.f, 1.f - z * z));
        float phi = 2 * Pi * u[1];
        return Tangent3f(r * std::cos(phi), r * std::sin(phi), z);
    }

    inline float UniformSpherePdf() { return Inv4Pi; }

    inline Tangent3f CosineSampleHemisphere(const Tuple2<float>& u) {
        Point2f d = ConcentricSampleDisk(u);
        float z = std::sqrt(std::max(0.f, 1 - d.x * d.x - d.y * d.y));
        return Tangent3f(d.x, d.y, z);
    }

    inline float CosineHemispherePdf(float cosTheta) { return cosTheta * InvPi; }

    inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return (f * f) / (f * f + g * g);
    }

} // namespace lightfold