#pragma once

#include <core/interaction.h>
#include <core/spectrum.h>
#include <math/sample.h>

namespace lightfold {

    class BSDF {
    public:
        // BSDF Interface
        virtual ~BSDF() {}
        BSDF(const SurfaceInteraction& si, bool isSingular) : ns(si.shading.n), ng(si.n),
            ss(Normalize(si.shading.dpdu)), ts(Cross(ns, ss)), isSingular(isSingular) {}

        virtual Spectrum distF(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        // Gives the value of the distribution function for the given pair of directions
        virtual Spectrum Sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& sample,
            float* pdf, bool isSingular = false) const = 0;
        // Samples an incident direction from the given outgoing direction and returns the value of the BxDF.
        virtual Spectrum reflectance(const Tangent3f& wo, int nSamples, const Point2f* samples) const;
        // Compute hemispherical-directional reflectance
        virtual Spectrum reflectance(int nSamples, const Point2f* samples1, const Point2f* samples2) const;
        // Compute hemispherical-hemispherical reflectance
        virtual float Pdf(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        // Gives the value of the PDF for the given pair of directions
        Tangent3f LocalToWorld(const Tangent3f& v) const {
            return Tangent3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                ss.y * v.x + ts.y * v.y + ns.y * v.z,
                ss.z * v.x + ts.z * v.y + ns.z * v.z);
        }

        // BSDF Public Data
        const bool isSingular;

    protected:
        // BSDF Protected Data
        const Normal3f ns, ng;
        const Tangent3f ss, ts;
    };

} // namespace lightfold