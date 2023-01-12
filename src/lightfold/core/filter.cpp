#include <core/filter.h>

namespace lightfold {

	float MitchellFilter::Evaluate(const Point2f& p) const {
		return Mitchell1D(p.x * invRadius.x) * Mitchell1D(p.y * invRadius.y);
	}

} // namespace lightfold