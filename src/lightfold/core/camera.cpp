#include <core/camera.h>

#include <core/stats.h>

namespace lightfold {

    float Camera::GenerateRayDifferential(const CameraSample& sample,
        RayDifferential* rd) const {
        float wt = GenerateRay(sample, rd);
        if (wt == 0) return 0;

        // Find camera ray after shifting a fraction of a pixel in the $x$ direction
        float wtx;
        for (float eps : { .05, -.05 }) {
            CameraSample sshift = sample;
            sshift.pFilm.x += eps;
            Ray rx;
            wtx = GenerateRay(sshift, &rx);
            rd->rxOrigin = rd->o + (rx.o - rd->o) / eps;
            rd->rxDirection = rd->d + (rx.d - rd->d) / eps;
            if (wtx != 0)
                break;
        }
        if (wtx == 0)
            return 0;

        // Find camera ray after shifting a fraction of a pixel in the $y$ direction
        float wty;
        for (float eps : { .05, -.05 }) {
            CameraSample sshift = sample;
            sshift.pFilm.y += eps;
            Ray ry;
            wty = GenerateRay(sshift, &ry);
            rd->ryOrigin = rd->o + (ry.o - rd->o) / eps;
            rd->ryDirection = rd->d + (ry.d - rd->d) / eps;
            if (wty != 0)
                break;
        }
        if (wty == 0)
            return 0;

        rd->hasDifferentials = true;
        return wt;
    }

    PerspectiveCamera::PerspectiveCamera(const AnimatedTransform& CameraToWorld,
        const Bounds2f& screenWindow, float shutterOpen, float shutterClose,
        float lensRadius, float focalDistance, float fov, Film* film, const Medium* medium)
        : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
            screenWindow, shutterOpen, shutterClose, lensRadius,
            focalDistance, film, medium) {
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

    float PerspectiveCamera::GenerateRayDifferential(const CameraSample& sample,
        RayDifferential* ray) const {
        ProfilePhase prof(Prof::GenerateCameraRay);
        // Compute raster and camera sample positions
        Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
        Point3f pCamera = RasterToCamera(pFilm);
        Tangent3f dir = Normalize(Tangent3f(pCamera.x, pCamera.y, pCamera.z));
        *ray = RayDifferential(Point3f(0, 0, 0), dir);
        // Modify ray for depth of field
        if (lensRadius > 0) {
            // Sample point on lens
            Point2f pLens = (Point2f)(lensRadius * (Vector2f)ConcentricSampleDisk(sample.pLens));

            // Compute point on plane of focus
            float ft = focalDistance / ray->d.z;
            Point3f pFocus = (*ray)(ft);

            // Update ray for effect of lens
            ray->o = Point3f(pLens.x, pLens.y, 0);
            ray->d = Normalize(pFocus - ray->o);
        }

        // Compute offset rays for _PerspectiveCamera_ ray differentials
        if (lensRadius > 0) {
            // Compute _PerspectiveCamera_ ray differentials accounting for lens

            // Sample point on lens
            Point2f pLens = (Point2f)(lensRadius * (Vector2f)ConcentricSampleDisk(sample.pLens));
            Tangent3f dx = Normalize(Tangent3f(pCamera.x, pCamera.y, pCamera.z)
                + Tangent3f(dxCamera.x, dxCamera.y, dxCamera.z));
            float ft = focalDistance / dx.z;
            Point3f pFocus = Point3f(0, 0, 0) + (ft * dx);
            ray->rxOrigin = Point3f(pLens.x, pLens.y, 0);
            ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

            Tangent3f dy = Normalize(Tangent3f(pCamera.x, pCamera.y, pCamera.z)
                + Tangent3f(dyCamera.x, dyCamera.y, dyCamera.z));
            ft = focalDistance / dy.z;
            pFocus = Point3f(0, 0, 0) + (ft * dy);
            ray->ryOrigin = Point3f(pLens.x, pLens.y, 0);
            ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
        }
        else {
            ray->rxOrigin = ray->ryOrigin = ray->o;
            ray->rxDirection = Normalize(Tangent3f(pCamera.x, pCamera.y, pCamera.z)
                + Tangent3f(dxCamera.x, dxCamera.y, dxCamera.z));
            ray->ryDirection = Normalize(Tangent3f(pCamera.x, pCamera.y, pCamera.z)
                + Tangent3f(dyCamera.x, dyCamera.y, dyCamera.z));
        }
        ray->time = std::lerp(sample.time, shutterOpen, shutterClose);
        ray->medium = medium;
        *ray = CameraToWorld(*ray);
        ray->hasDifferentials = true;
        return 1;
    }

} // namespace lightfold