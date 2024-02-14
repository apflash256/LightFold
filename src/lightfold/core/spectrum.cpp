#include <core/spectrum.h>

namespace lightfold {

#if SPECTRUM == GERYSCALE
	RGB Spectrum::GetRGB() const {
		return RGB(val, val, val);
	}
#endif

} // namespace lightfold