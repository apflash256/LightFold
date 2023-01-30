#pragma once

#include <core/interaction.h>
#include <core/spectrum.h>



namespace lightfold {

    class BSSRDF {
    public:
        // BSSRDF Public Methods
        BSSRDF(const SurfaceInteraction& po, float eta) : po(po), eta(eta) {}
        virtual ~BSSRDF() {}

        // BSSRDF Interface
        virtual Spectrum S(const SurfaceInteraction& pi, const Tangent3f& wi) = 0;
        virtual Spectrum Sample_S(const Scene& scene, float u1, const Point2f& u2,
            SurfaceInteraction* si, float* pdf) const = 0;

    protected:
        // BSSRDF Protected Data
        const SurfaceInteraction& po;
        float eta;
    };

} // namespace lightfold