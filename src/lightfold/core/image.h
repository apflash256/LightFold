#pragma once

#include <memory>

namespace lightfold {

	int WriteEXR(std::unique_ptr<float[]> rgb, int width, int height, const char* outfilename);

} // namespace lightfold