#pragma once
#include <util/vector.h>
#include <core/film.h>

namespace lightfold {

	class CameraSample {
	public:
		// CameraSample Public Methods
		CameraSample() : filmSample(Point2<Float>()), lensSample(Point2<Float>()), time(0) {}
		CameraSample(Point2<Float> filmSample, Point2<Float> lensSample, Float time) :
			filmSample(filmSample), lensSample(lensSample), time(time) {}

		// CameraSample Public Members
		Point2<Float> filmSample; // (0,0) to (width, height)
		Point2<Float> lensSample; // (-1,-1) to (1,1)
		Float time;
	};

	class Camera {
	public:
		// Camera Public Methods
		Camera(Point origin, Transformation transformation, Film film) :
			origin(origin), transformation(transformation), film(film) {}

		virtual Ray GenerateRay(CameraSample cameraSample) const {} = 0;

	protected:
		// Camera Protected Members
		Point origin;
		Transformation transformation;
		Film film;
	};

} // namespace lightfold