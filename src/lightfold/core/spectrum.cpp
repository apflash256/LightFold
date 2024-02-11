#include <core/spectrum.h>

namespace lightfold {

#if SPECTRUM == GERYSCALE
	RGB Spectrum::getRGB() const {
		struct RGB rgb = { val, val, val };
		return rgb;
	}
#endif

} // namespace lightfold