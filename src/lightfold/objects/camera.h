#pragma once

#include <core/film.h>
#include <math/sample.h>
#include <math/transform.h>

namespace lightfold {

	struct CameraSample {
		Point2f pFilm;
		Point2f pLens;
	};

	class Camera {
	public:
		// Camera Public Methods
		Camera(const Transform& CameraToWorld, Film* film) : CameraToWorld(CameraToWorld), film(film) {}
		virtual float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;

		// Camera Public Data
		Transform CameraToWorld;
		Film* film;
	};

	class ProjectiveCamera : public Camera {
	public:
		// ProjectiveCamera Public Methods
		ProjectiveCamera(const Transform& CameraToWorld,
			const Transform& CameraToScreen, const Bounds2f& screenWindow, float lensr, float focald, Film* film)
			: Camera(CameraToWorld, film), CameraToScreen(CameraToScreen) {
			lensRadius = lensr;
			focalDistance = focald;
			ScreenToRaster = Scale((float)film->fullResolution.x, (float)film->fullResolution.y, 1) *
				Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x), 1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
				Translate(TanVector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
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
		PerspectiveCamera(const Transform& CameraToWorld, const Bounds2f& screenWindow,
			float lensRadius, float focalDistance, float fov, Film* film);
		float GenerateRay(const CameraSample& sample, Ray*) const;

	private:
		// PerspectiveCamera Private Data
		TanVector3f dxCamera, dyCamera;
		float A;
	};

} // namespace lightfold