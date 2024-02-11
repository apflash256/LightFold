#pragma once
#include "core/image.h"

#include <memory>

using namespace lightfold;

constexpr int width = 1280;
constexpr int height = 720;

int main(void) {
	char filename[] = "image.exr";
	std::unique_ptr<float[]> img(new float[width * height * 3]);
	WriteEXR(img, width, height, filename);
	return 0;
}