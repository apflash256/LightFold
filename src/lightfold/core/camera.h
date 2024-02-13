#pragma once
#include <util/vector.h>

namespace lightfold {

	class CameraSample {
	public:
		// CameraSample Public Methods
		CameraSample() : filmSample(Point2<Float>()), lensSample(Point2<Float>()), time(0) {}
		CameraSample(Point2<Float> filmSample, Point2<Float> lensSample, Float time) :
			filmSample(filmSample), lensSample(lensSample), time(time) {}

		// CameraSample Public Members
		Point2<Float> filmSample, lensSample;
		Float time;
	};

} // namespace lightfold