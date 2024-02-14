#include <camera/perspective.h>
#include <algorithm>

namespace lightfold {

	Ray PerspectiveCamera::GenerateRay(CameraSample cameraSample) const {
		Vector startPoint(aperture * cameraSample.lensSample.x,
			aperture * cameraSample.lensSample.y, 0);
		int maxForm = std::max(film.width, film.height);
		Vector endPoint = (tanTheta / maxForm) *
			Vector(film.width / 2 - cameraSample.filmSample.x,
				cameraSample.filmSample.y - film.height / 2, 1);
		return Ray(Transport(origin, startPoint), Transform(endPoint - startPoint, transformation));
	}

} // namespace lightfold