#include <core/scene.h>
#include <aggregate/bvh.h>
#include <utils/uvsphere.h>

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

	Point3f pos(1, 2, 1);
	Point3f look(0, 0, 0);
	Tangent3f up(0, 0, 1);
	Transform c2w = Inverse(LookAt(pos, look, up));
	Bounds2f screenWindow = { {-((float)WIDTH) / HEIGHT, -1.f},{((float)WIDTH) / HEIGHT, 1.f} };
	float lensRadius = 0.f, focalDistance = 0.04f, fov = 1.0f;
	PerspectiveCamera myCam(c2w, screenWindow, lensRadius, focalDistance, fov, &myFilm);

	Bounds2i sampleBounds = myCam.film->GetSampleBounds();
	Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize, (sampleExtent.y + tileSize - 1) / tileSize);

	auto prims = UVSphere(12, 24);
	std::shared_ptr<Primitive> bvh = std::make_shared<BVHAccel>(prims);

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
					SurfaceInteraction isect;
					float depth = 0;
					if (bvh->Intersect(ray, &isect)) {
						depth = 1.f / (Length(isect.p - pos) - 1.4);
						//std::cout << depth << std::endl;
					}
					filmTile->AddSample(cameraSample.pFilm, Spectrum(depth), rayWeight);
				} while (tileSampler->StartNextSample());
			}
			myCam.film->MergeFilmTile(std::move(filmTile));
		}, nTiles);
	myCam.film->WriteImage();
	ParallelCleanup();

	return 0;
}