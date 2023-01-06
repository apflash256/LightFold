#include <iostream>
#include <imageio.h>
#include <math/transform.h>
#include <objects/camera.h>
#include <objects/shape.h>

#define WIDTH 3840
#define HEIGHT 2160

using namespace lightfold;
float testimg[3 * WIDTH * HEIGHT];

int main(void) {
	Transform idt = Scale(1.f, 1.f, 1.f);
	Bounds2f screenWindow = { {-((float)WIDTH) / HEIGHT, -1.f},{((float)WIDTH) / HEIGHT, 1.f} };
	float lensRadius = 0.f, focalDistance = 0.04f, fov = 1.0f;
	Film uhd = { {WIDTH, HEIGHT}, 1.f };
	PerspectiveCamera myCam(idt, screenWindow, lensRadius, focalDistance, fov, &uhd);

	Transform otw = Scale(1.f, 1.f, 1.f);
	Point3f p[3] = { { -1.0f, -1.0f, 2.f } ,{ 1.0f, -1.0f, 4.f } ,{ 1.f, 1.0f, 6.f } };
	TanVector3f* s = nullptr;
	CotVector3f* n = nullptr;
	Point2f* uv = nullptr;
	int vInds[3] = { 0,1,2 };
	std::shared_ptr<TriangleMesh> myMesh = std::make_shared<TriangleMesh>(otw, 1, vInds, 3, p, s, n, uv);

	Transform wto = Scale(1.f, 1.f, 1.f);
	Triangle myTri(&otw, &wto, false, myMesh, 0);

	float f, depth;
	CameraSample camSample;
	Ray ray;
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < HEIGHT; j++) {
			camSample = { {(float)i,(float)j},{0.f,0.f} };
			myCam.GenerateRay(camSample, &ray);
			depth = myTri.Intersect(ray, &f) ? 1/f : 0.f;
			testimg[3 * i + 3 * WIDTH * j] = depth;
			testimg[3 * i + 3 * WIDTH * j + 1] = depth;
			testimg[3 * i + 3 * WIDTH * j + 2] = depth;
		}
	}
	int a = WriteEXR(testimg, WIDTH, HEIGHT, "testimg.exr");
	return 0;
}