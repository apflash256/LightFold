#include <core/spectrum.h>

namespace lightfold {

#if SPECTRUM == GERYSCALE
	RGB Spectrum::getRGB() const {
		return RGB(val, val, val);
	}
#endif

} // namespace lightfold