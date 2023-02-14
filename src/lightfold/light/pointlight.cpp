#include <light/pointlight.h>

#include <core/stats.h>
#include <math/sample.h>

namespace lightfold {

    Spectrum PointLight::Sample_Li(const Interaction& ref, const Point2f& u,
        Tangent3f* wi, float* pdf,
        VisibilityTester* vis) const {
        ProfilePhase _(Prof::LightSample);
        *wi = Normalize(pLight - ref.p);
        *pdf = 1.f;
        *vis = VisibilityTester(ref, Interaction(pLight, ref.time, mediumInterface));
        return I / DistanceSquared(pLight, ref.p);
    }

    Spectrum PointLight::Power() const { return 4 * Pi * I; }

    float PointLight::Pdf_Li(const Interaction&, const Tangent3f&) const {
        return 0;
    }

    Spectrum PointLight::Sample_Le(const Point2f& u1, const Point2f& u2, float time,
        Ray* ray, Normal3f* nLight, float* pdfPos,
        float* pdfDir) const {
        ProfilePhase _(Prof::LightSample);
        *ray = Ray(pLight, UniformSampleSphere(u1), Infinity, time,
            mediumInterface.inside);
        *nLight = (Normal3f)ray->d;
        *pdfPos = 1;
        *pdfDir = UniformSpherePdf();
        return I;
    }

    void PointLight::Pdf_Le(const Ray&, const Normal3f&, float* pdfPos,
        float* pdfDir) const {
        ProfilePhase _(Prof::LightPdf);
        *pdfPos = 0;
        *pdfDir = UniformSpherePdf();
    }

} // namespace lightfold