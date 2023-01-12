#include <iostream>
#include <core/imageio.h>
#include <core/spectra.h>
#include <core/filter.h>
#include <core/film.h>
#include <math/transform.h>
#include <objects/camera.h>
#include <objects/shape.h>
#include <sample/sampler.h>

#define WIDTH 3840
#define HEIGHT 2160

using namespace lightfold;
float testimg[3 * WIDTH * HEIGHT];

int main(void) {
	Point2i uhd(WIDTH, HEIGHT);
	Bounds2f cropWindow({ 0.f,0.f }, { 1.f, 1.f });
	std::unique_ptr<Filter> myFilter(new MitchellFilter({ 2.f,2.f }, 1.f / 3.f, 1.f / 3.f));
	char filename[] = "testimg.exr";
	Film myFilm(uhd, cropWindow, std::move(myFilter), 1.f, filename, 1.f);
	Transform idt = Scale(1.f, 1.f, 1.f);
	Bounds2f screenWindow = { {-((float)WIDTH) / HEIGHT, -1.f},{((float)WIDTH) / HEIGHT, 1.f} };
	float lensRadius = 0.f, focalDistance = 0.04f, fov = 1.0f;
	PerspectiveCamera myCam(idt, screenWindow, lensRadius, focalDistance, fov, &myFilm);

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
	Point2f pFilm;
	CameraSample camSample;
	Ray ray;
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < HEIGHT; j++) {
			pFilm.x = (float)i;
			pFilm.y = (float)j;
			camSample = { pFilm, {0.f,0.f} };
			myCam.GenerateRay(camSample, &ray);
			depth = myTri.Intersect(ray, &f) ? 1/f : 0.f;
			myFilm.AddSplat(pFilm, BW(depth));
		}
	}

	myFilm.WriteImage();

	return 0;
}