#pragma once

#include <core/rgb.h>

#include <memory>

namespace lightfold {

	int WriteEXR(std::unique_ptr<RGB[]> rgb, int width, int height, const char* outfilename);

} // namespace lightfold