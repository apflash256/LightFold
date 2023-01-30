#pragma once

#include <core/sampler.h>

namespace lightfold {

    class MediumInteraction;
    class Ray;

    class PhaseFunction {
    public:
        // PhaseFunction Interface
        virtual ~PhaseFunction() { }
        virtual float p(const Tangent3f& wo, const Tangent3f& wi) const = 0;
        virtual float Sample_p(const Tangent3f& wo, Tangent3f* wi, const Point2f& u) const = 0;
    };

    class Medium {
    public:
        // Medium Interface
        virtual ~Medium() { }
        virtual Spectrum Transmittance(const Ray& ray, Sampler& sampler) const = 0;
        virtual Spectrum Sample(const Ray& ray, Sampler& sampler, MediumInteraction* mi) const = 0;
    };

    struct MediumInterface {
        // MediumInterface Public Methods
        MediumInterface() : inside(nullptr), outside(nullptr) { }
        MediumInterface(const Medium * medium) : inside(medium), outside(medium) { }
        MediumInterface(const Medium* inside, const Medium* outside) : inside(inside), outside(outside) { }
        bool IsMediumTransition() const { return inside != outside; }

        // MediumInterface Public Data
        const Medium* inside, * outside;
    };

} // namespace lightfold