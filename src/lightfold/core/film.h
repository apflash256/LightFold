#pragma once
#include <core/settings.h>
#include <core/rgb.h>

namespace lightfold {

	class Film {
	public:
		// Film Public Methods
		RGB GetPixel() const;

		// Film Public Members
		int width, height;
	};

} // namespace lightfold