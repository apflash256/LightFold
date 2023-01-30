#pragma once

namespace lightfold {

	class SurfaceInteraction;

	template <typename T> class Texture {
	public:
		// SurfaceTexture Methods
		virtual T Evaluate(const SurfaceInteraction&) const = 0;
		virtual ~Texture() { }
	};

} // namespace lightfold