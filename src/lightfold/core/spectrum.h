#pragma once
#include <core/settings.h>
#include <core/rgb.h>

namespace lightfold {

	class Spectrum {
	public:
		RGB getRGB() const;

	private:
#if SPECTRUM == GREYSCALE
		Float val;
#endif
	};

} // namespace lightfold