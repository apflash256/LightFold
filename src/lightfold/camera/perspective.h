#pragma once
#include <core/camera.h>

namespace lightfold {

	class PerspectiveCamera : public Camera {
	public:
		// PerspectiveCamera Public Methods
		PerspectiveCamera(Point origin, Transformation transformation, Film film,
			Float focalDistance, Float aperture, Float tanTheta) :
			origin(origin), transformation(transformation), film(film),
			focalDistance(focalDistance), aperture(aperture), tanTheta(tanTheta) {}

	private:
		// PerspectiveCamera Private Members
		Float focalDistance, aperture, tanTheta;
	};

} // namespace lightfold
