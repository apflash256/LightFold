#pragma once

#include <core/lobe.h>
#include <core/interaction.h>

namespace lightfold {

    class BSDF {
    public:
        // BSDF Public Methods
        BSDF(const SurfaceInteraction& si)
            : ns(si.shading.n), ng(si.n), ss(Normalize(si.shading.dpdu)), ts(Cross(ns, ss)) {}
        ~BSDF() { }
        void Add(std::shared_ptr<Lobe> b) {
            lobes.push_back(b);
        }
        int NumComponents(LobeType flags = LOBE_ALL) const {
            int num = 0;
            for (int i = 0; i < size(lobes); ++i)
                if (lobes[i]->matchesFlags(flags)) ++num;
            return num;
        }
        Tangent3f worldToLocal(const Tangent3f& v) const {
            return Tangent3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
        }
        Tangent3f localToWorld(const Tangent3f& v) const {
            return Tangent3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                ss.y * v.x + ts.y * v.y + ns.y * v.z,
                ss.z * v.x + ts.z * v.y + ns.z * v.z);
        }
        Spectrum dist_F(const Tangent3f& woW, const Tangent3f& wiW,
            LobeType flags = LOBE_ALL) const;
        Spectrum reflectance(int nSamples, const Point2f* samples1, const Point2f* samples2,
            LobeType flags = LOBE_ALL) const;
        Spectrum reflectance(const Tangent3f& wo, int nSamples, const Point2f* samples,
            LobeType flags = LOBE_ALL) const;
        Spectrum sample_F(const Tangent3f& wo, Tangent3f* wi, const Point2f& u,
            float* pdf, LobeType type = LOBE_ALL, LobeType* sampledType = nullptr) const;
        float PDF(const Tangent3f& wo, const Tangent3f& wi, LobeType flags = LOBE_ALL) const;

    private:
        // BSDF Private Data
        const Normal3f ns, ng;
        const Tangent3f ss, ts;
        std::vector<std::shared_ptr<Lobe>> lobes;
        friend class MixMaterial;
    };

} // namespace lightfold