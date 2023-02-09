#pragma once

#include <core/film.h>
#include <math/sample.h>
#include <math/transform.h>

namespace lightfold {

	struct CameraSample {
		Point2f pFilm;
		Point2f pLens;
		float time;
	};

	class Camera {
	public:
		// Camera Public Methods
		Camera(const AnimatedTransform& CameraToWorld, float shutterOpen, float shutterClose,
			Film* film, const Medium* medium)
			: CameraToWorld(CameraToWorld), shutterOpen(shutterOpen),
			shutterClose(shutterClose), film(film), medium(medium) {}
		virtual float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;
		virtual float GenerateRayDifferential(const CameraSample& sample,
			RayDifferential* rd) const;

		// Camera Public Data
		AnimatedTransform CameraToWorld;
		Film* film;
		const float shutterOpen, shutterClose;
		const Medium* medium;
	};

	class ProjectiveCamera : public Camera {
	public:
		// ProjectiveCamera Public Methods
		ProjectiveCamera(const AnimatedTransform& CameraToWorld,
			const Transform& CameraToScreen, const Bounds2f& screenWindow, float shutterOpen,
			float shutterClose, float lensr, float focald, Film* film, const Medium* medium)
			: Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
			CameraToScreen(CameraToScreen) {
			lensRadius = lensr;
			focalDistance = focald;
			ScreenToRaster = Scale((float)film->fullResolution.x, (float)film->fullResolution.y, 1) *
				Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x), 1 / (screenWindow.pMax.y - screenWindow.pMin.y), 1) *
				Translate(Tangent3f(-screenWindow.pMin.x, -screenWindow.pMin.y, 0));
			RasterToScreen = Inverse(ScreenToRaster);
			RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
		}

	protected:
		// ProjectiveCamera Protected Data
		Transform CameraToScreen, RasterToCamera;
		Transform ScreenToRaster, RasterToScreen;
		float lensRadius, focalDistance;
	};

	class PerspectiveCamera : public ProjectiveCamera {
	public:
		// PerspectiveCamera Public Methods
		PerspectiveCamera(const AnimatedTransform& CameraToWorld, const Bounds2f& screenWindow,
			float shutterOpen, float shutterClose, float lensRadius, float focalDistance,
			float fov, Film* film, const Medium* medium);
		float GenerateRay(const CameraSample& sample, Ray*) const;
		float GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const;

	private:
		// PerspectiveCamera Private Data
		Tangent3f dxCamera, dyCamera;
		float A;
	};

} // namespace lightfold