#include <core/camera.h>

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
        A = std::abs((pMax.x / pMax.z - pMin.x / pMin.z) * (pMax.y / pMax.z - pMin.y / pMin.z));
    }

    float PerspectiveCamera::GenerateRay(const CameraSample& sample,
        Ray* ray) const {
        Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
        Point3f pCamera = RasterToCamera(pFilm);
        *ray = Ray(Point3f(0, 0, 0), Normalize(pCamera - Point3f(0, 0, 0)));
        if (lensRadius > 0) {
            Point2f pLens = (Point2f)(lensRadius * (Vector2f)ConcentricSampleDisk(sample.pLens));
            float ft = focalDistance / ray->d.z;
            Point3f pFocus = (*ray)(ft);
            ray->o = Point3f(pLens.x, pLens.y, 0);
            ray->d = Normalize(pFocus - ray->o);
        }
        *ray = CameraToWorld(*ray);
        return 1.f;
    }

} // namespace lightfold