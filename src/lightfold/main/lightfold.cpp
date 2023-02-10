#include <core/scene.h>
#include <aggregate/bvh.h>
#include <light/pointlight.h>
#include <material/matte.h>
#include <texture/constanttexture.h>
#include <utils/uvsphere.h>
#include <integrator/directlighting.h>

constexpr auto WIDTH = 3840;
constexpr auto HEIGHT = 2160;
constexpr auto SPP = 10;

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
	AnimatedTransform c2wa(c2w, 0, c2w, 1);
	Bounds2f screenWindow = { {-((float)WIDTH) / HEIGHT, -1.f},{((float)WIDTH) / HEIGHT, 1.f} };
	float lensRadius = 0.f, focalDistance = 0.04f, fov = 1.0f;
	std::shared_ptr<Camera> myCam = std::make_shared<PerspectiveCamera>(c2wa, screenWindow,
		0, 0, lensRadius, focalDistance, fov, &myFilm, nullptr);

	Bounds2i sampleBounds = myCam->film->GetSampleBounds();
	Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize, (sampleExtent.y + tileSize - 1) / tileSize);

	std::shared_ptr<Texture<Spectrum>> color =
		std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.6, 0.6, 0.6));
	std::shared_ptr<Texture<float>> roughness =
		std::make_shared<ConstantTexture<float>>(0.4f);
	std::shared_ptr<Material> myMaterial =
		std::make_shared<MatteMaterial>(color, roughness);
	auto prims = UVSphere(120, 240, myMaterial, nullptr, MediumInterface());
	std::shared_ptr<Primitive> bvh = std::make_shared<BVHAccel>(prims);

	Transform ltw = Translate(Tangent3f(4, 1, 1));
	std::shared_ptr<Light> myLight = std::make_shared<PointLight>(ltw, nullptr, Spectrum(10, 10, 10));
	std::vector<std::shared_ptr<Light>> myLights = { myLight };

	Scene myScene(bvh, myLights);

	std::shared_ptr<Sampler> mySampler = std::make_shared<HaltonSampler>(SPP, sampleBounds);

	DirectLightingIntegrator myIntegrator(LightStrategy::UniformSampleAll, 1, myCam, mySampler, Bounds2i({ 0, 0 }, uhd));
	myIntegrator.Render(myScene);

	return 0;
}