#pragma once

#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

int WriteEXR(const float* rgb, int width, int height, const char* outfilename);