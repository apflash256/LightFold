#include <core/bsdf.h>

#include <core/stats.h>

namespace lightfold {

    Spectrum BSDF::dist_F(const Tangent3f& woW, const Tangent3f& wiW,
        LobeType flags) const {
        ProfilePhase pp(Prof::BSDFEvaluation);
        Tangent3f wi = worldToLocal(wiW), wo = worldToLocal(woW);
        if (wo.z == 0) return 0.;
        bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
        Spectrum f(0.f);
        for (int i = 0; i < size(lobes); ++i)
            if (lobes[i]->matchesFlags(flags) &&
                ((reflect && (lobes[i]->type & LOBE_REFLECTION)) ||
                    (!reflect && (lobes[i]->type & LOBE_TRANSMISSION))))
                f += lobes[i]->dist_F(wo, wi);
        return f;
    }

    Spectrum BSDF::reflectance(int nSamples, const Point2f* samples1,
        const Point2f* samples2, LobeType flags) const {
        Spectrum ret(0.f);
        for (int i = 0; i < size(lobes); ++i)
            if (lobes[i]->matchesFlags(flags))
                ret += lobes[i]->reflectance(nSamples, samples1, samples2);
        return ret;
    }

    Spectrum BSDF::reflectance(const Tangent3f& woWorld, int nSamples, const Point2f* samples,
        LobeType flags) const {
        Tangent3f wo = worldToLocal(woWorld);
        Spectrum ret(0.f);
        for (int i = 0; i < size(lobes); ++i)
            if (lobes[i]->matchesFlags(flags))
                ret += lobes[i]->reflectance(wo, nSamples, samples);
        return ret;
    }

    Spectrum BSDF::sample_F(const Tangent3f& woWorld, Tangent3f* wiWorld, const Point2f& u,
        float* pdf, LobeType type, LobeType* sampledType) const {
        ProfilePhase pp(Prof::BSDFSampling);
        // Choose which lobe to sample
        int matchingComps = NumComponents(type);
        if (matchingComps == 0) {
            *pdf = 0;
            if (sampledType) *sampledType = LobeType(0);
            return Spectrum(0);
        }
        int comp = std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

        // Get lobe pointer for chosen component
        std::shared_ptr<Lobe> lobe = nullptr;
        int count = comp;
        for (int i = 0; i < size(lobes); ++i)
            if (lobes[i]->matchesFlags(type) && count-- == 0) {
                lobe = lobes[i];
                break;
            }

        // Remap lobe sample _u_ to $[0,1)^2$
        Point2f uRemapped(std::min(u[0] * matchingComps - comp, OneMinusEpsilon), u[1]);

        // Sample chosen lobe
        Tangent3f wi, wo = worldToLocal(woWorld);
        if (wo.z == 0) return 0.;
        *pdf = 0;
        if (sampledType) *sampledType = lobe->type;
        Spectrum f = lobe->sample_F(wo, &wi, uRemapped, pdf, sampledType);
        if (*pdf == 0) {
            if (sampledType) *sampledType = LobeType(0);
            return 0;
        }
        *wiWorld = localToWorld(wi);

        // Compute overall PDF with all matching lobes
        if (!(lobe->type & LOBE_SPECULAR) && matchingComps > 1)
            for (int i = 0; i < size(lobes); ++i)
                if (lobes[i] != lobe && lobes[i]->matchesFlags(type))
                    *pdf += lobes[i]->PDF(wo, wi);
        if (matchingComps > 1) *pdf /= matchingComps;

        // Compute value of lobe for sampled direction
        if (!(lobe->type & LOBE_SPECULAR)) {
            bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
            f = 0.;
            for (int i = 0; i < size(lobes); ++i)
                if (lobes[i]->matchesFlags(type) &&
                    ((reflect && (lobes[i]->type & LOBE_REFLECTION)) ||
                        (!reflect && (lobes[i]->type & LOBE_TRANSMISSION))))
                    f += lobes[i]->dist_F(wo, wi);
        }
        return f;
    }

    float BSDF::PDF(const Tangent3f& woWorld, const Tangent3f& wiWorld, LobeType flags) const {
        ProfilePhase pp(Prof::BSDFPdf);
        if (size(lobes) == 0.f) return 0.f;
        Tangent3f wo = worldToLocal(woWorld), wi = worldToLocal(wiWorld);
        if (wo.z == 0) return 0.;
        float pdf = 0.f;
        int matchingComps = 0;
        for (int i = 0; i < size(lobes); ++i)
            if (lobes[i]->matchesFlags(flags)) {
                ++matchingComps;
                pdf += lobes[i]->PDF(wo, wi);
            }
        float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
        return v;
    }

} // namespace lightfold