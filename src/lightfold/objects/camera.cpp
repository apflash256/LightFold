#include <objects/camera.h>

namespace lightfold {

    PerspectiveCamera::PerspectiveCamera(
        const Transform& CameraToWorld, const Bounds2f& screenWindow,
        float lensRadius, float focalDistance, float fov, Film* film)
        : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
            screenWindow, lensRadius, focalDistance, film) {
        dxCamera = (RasterToCamera(Point3f(1, 0, 0)) - RasterToCamera(Point3f(0, 0, 0)));
        dyCamera = (RasterToCamera(Point3f(0, 1, 0)) - RasterToCamera(Point3f(0, 0, 0)));
        Point2i res = film->fullResolution;
        Point3f pMin = RasterToCamera(Point3f(0, 0, 0));
        Point3f pMax = RasterToCamera(Point3f((float)res.x, (float)res.y, 0));
        pMin /= pMin.z;
        pMax /= pMax.z;
        A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));
    }

    float PerspectiveCamera::GenerateRay(const CameraSample& sample,
        Ray* ray) const {
        Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
        Point3f pCamera = RasterToCamera(pFilm);
        *ray = Ray(Point3f(0, 0, 0), Normalize((TanVector3f)(GeoVector3f)pCamera));
        if (lensRadius > 0) {
            Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
            float ft = focalDistance / ray->d.z;
            Point3f pFocus = (*ray)(ft);
            ray->o = Point3f(pLens.x, pLens.y, 0);
            ray->d = Normalize(pFocus - ray->o);
        }
        *ray = CameraToWorld(*ray);
        return 1;
    }

} // namespace lightfold