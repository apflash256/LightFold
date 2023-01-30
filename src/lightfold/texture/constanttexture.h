#pragma once

#include <core/texture.h>

namespace lightfold {

	template <typename T> class ConstantTexture : public Texture<T> {
	public:
		// ConstantTexture Public Methods
		T Evaluate(const SurfaceInteraction&) const { return value; }

		// ConstantTexture Public Data
		T value;
	};

} // namespace lightfold