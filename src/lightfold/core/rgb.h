#include <core/settings.h>

namespace lightfold {

	class RGB {
	public:
		// RGB Public Methods
		RGB() : r(0), g(0), b(0) {}
		RGB(Float r, Float g, Float b) : r(r), g(g), b(b) {}

		// RGB Public Members
		Float r, g, b;
	};

} // namespace lightfold