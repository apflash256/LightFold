#pragma once
#include "core/image.h"

#include <memory>

using namespace lightfold;

constexpr int width = 1280;
constexpr int height = 720;

int main(void) {
	char filename[] = "image.exr";
	std::unique_ptr<RGB[]> img(new RGB[width * height]);
	WriteEXR(std::move(img), width, height, filename);
	return 0;
}