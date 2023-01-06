#include <iostream>
#include <imageio.h>
#include <math/transform.h>
#include <objects/camera.h>
#include <objects/shape.h>

using namespace lightfold;

int main(void) {
	Transform idt = Scale(1.f, 1.f, 1.f);
	Bounds2f screenWindow = { {-0.384f, -0.216f},{0.384f, 0.216f} };
	float lensRadius = 0.02f, focalDistance = 0.25f, fov = 1.6f;
	Film uhd = { {3840, 2160}, 1.f };
	PerspectiveCamera A(idt,screenWindow,lensRadius,focalDistance,fov,&uhd);
	std::cout << "Hello World!";
	float testimg[30000];
	for (int i = 0; i < 30000; i++) {
		testimg[i] = 1.0f * i;
	}
	int a = WriteEXR(testimg, 100, 100, "testimg.exr");
	return a;
}