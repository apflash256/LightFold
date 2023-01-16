#include <core/film.h>
#include <objects/camera.h>
#include <objects/shape.h>
#include <sample/sampler.h>

constexpr auto WIDTH = 3840;
constexpr auto HEIGHT = 2160;
constexpr auto SPP = 16;

using namespace lightfold;

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

	Bounds2i sampleBounds = myCam.film->GetSampleBounds();
	GeoVector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize, (sampleExtent.y + tileSize - 1) / tileSize);

	Transform otw = Scale(1.f, 1.f, 1.f);
	Point3f p[3] = { { -1.0f, -1.0f, 2.f } ,{ 1.0f, -1.0f, 4.f } ,{ 1.f, 1.0f, 6.f } };
	TanVector3f* s = nullptr;
	CotVector3f* n = nullptr;
	Point2f* uv = nullptr;
	int vInds[3] = { 0,1,2 };
	std::shared_ptr<TriangleMesh> myMesh = std::make_shared<TriangleMesh>(otw, 1, vInds, 3, p, s, n, uv);
	Transform wto = Scale(1.f, 1.f, 1.f);
	Triangle myTri(&otw, &wto, false, myMesh, 0);

	HaltonSampler mySampler(SPP, sampleBounds);
	ParallelInit();
	ParallelFor2D(
		[&](Point2i tile) {
			int seed = tile.y * nTiles.x + tile.x;
			std::unique_ptr<Sampler> tileSampler = mySampler.Clone(seed);

			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
			int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
			int y0 = sampleBounds.pMin.y + tile.y * tileSize;
			int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
			Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

			std::unique_ptr<FilmTile> filmTile = myCam.film->GetFilmTile(tileBounds);
			for (Point2i pixel : tileBounds) {
				tileSampler->StartPixel(pixel);
				do {
					CameraSample cameraSample = tileSampler->GetCameraSample(pixel);
					Ray ray;
					float rayWeight = myCam.GenerateRay(cameraSample, &ray);
					float f;
					float depth = myTri.Intersect(ray, &f) ? 1 / f : 0.f;
					Spectrum L(depth);
					// L is the data!
					filmTile->AddSample(cameraSample.pFilm, L, rayWeight);
				} while (tileSampler->StartNextSample());
			}
			myCam.film->MergeFilmTile(std::move(filmTile));
		}, nTiles);
	myCam.film->WriteImage();
	ParallelCleanup();
	return 0;
}