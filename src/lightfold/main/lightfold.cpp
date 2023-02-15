#include <core/scene.h>
#include <aggregate/bvh.h>
#include <light/pointlight.h>
#include <material/matte.h>
#include <shape/triangle.h>
#include <texture/constanttexture.h>
#include <utils/uvsphere.h>
#include <integrator/pathtracing.h>

#include <memory>

constexpr auto WIDTH = 1920;
constexpr auto HEIGHT = 1080;
constexpr auto SPP = 256;

using namespace lightfold;

int main(void) {
	Point2i uhd(WIDTH, HEIGHT);
	Bounds2f cropWindow({ 0.f, 0.f }, { 1.f, 1.f });
	std::unique_ptr<Filter> myFilter(new MitchellFilter({ 2.f,2.f }, 1.f / 3.f, 1.f / 3.f));
	char filename[] = "testimg.exr";
	Film myFilm(uhd, cropWindow, std::move(myFilter), 1.f, filename, 1.f);

	Point3f pos(0, -5, 0);
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

	Transform sphereotw = Translate({ 0, 0, -1.5 }) * Scale(1.5, 1.5, 1.5);
	Transform spherewto = Inverse(sphereotw);
	std::shared_ptr<Texture<Spectrum>> color1 =
		std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.8, 0.8, 0.8));
	std::shared_ptr<Texture<Spectrum>> color2 =
		std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.8, 0.1, 0.1));
	std::shared_ptr<Texture<Spectrum>> color3 =
		std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.1, 0.8, 0.1));
	std::shared_ptr<Texture<float>> roughness =
		std::make_shared<ConstantTexture<float>>(0.4f);
	std::shared_ptr<Material> myMaterial1 =
		std::make_shared<MatteMaterial>(color1, roughness);
	std::shared_ptr<Material> myMaterial2 =
		std::make_shared<MatteMaterial>(color2, roughness);
	std::shared_ptr<Material> myMaterial3 =
		std::make_shared<MatteMaterial>(color3, roughness);

	auto prims = UVSphere(&sphereotw, &spherewto, 120, 240, myMaterial1, nullptr, MediumInterface());

	Point3f boxPoints1[8] = { {-3, 3, 3}, {-3, -3, 3}, {3, 3, 3}, {3, -3, 3},
		{-3, 3, -3}, {-3, -3, -3}, {3, 3, -3}, {3, -3, -3} };
	Point3f boxPoints2[4] = { {-3, -3, 3}, {-3, 3, 3}, {-3, 3, -3}, {-3, -3, -3} };
	Point3f boxPoints3[4] = { {3, -3, 3}, {3, 3, 3}, {3, 3, -3}, {3, -3, -3} };
	int boxVInds1[18] = { 0, 2, 1, 1, 2, 3, 0, 4, 2, 2, 4, 6, 4, 5, 6, 5, 7, 6 };
	int boxVInds2[6] = { 0, 2, 1, 2, 0, 3 };
	int boxVInds3[6] = { 0, 1, 2, 2, 3, 0 };
	Transform iden = Scale(1, 1, 1);
	std::vector<std::shared_ptr<Shape>> MyBox1 = CreateTriangleMesh(&iden, &iden, false, 6,
		boxVInds1, 8, boxPoints1, nullptr, nullptr, nullptr);
	std::vector<std::shared_ptr<Shape>> MyBox2 = CreateTriangleMesh(&iden, &iden, false, 2,
		boxVInds2, 4, boxPoints2, nullptr, nullptr, nullptr);
	std::vector<std::shared_ptr<Shape>> MyBox3 = CreateTriangleMesh(&iden, &iden, false, 2,
		boxVInds3, 4, boxPoints3, nullptr, nullptr, nullptr);

	for (std::shared_ptr<Shape> shape : MyBox1) {
		prims.push_back(std::make_shared<GeometricPrimitive>(shape, myMaterial1, nullptr,
			MediumInterface()));
	}
	for (std::shared_ptr<Shape> shape : MyBox2) {
		prims.push_back(std::make_shared<GeometricPrimitive>(shape, myMaterial2, nullptr,
			MediumInterface()));
	}
	for (std::shared_ptr<Shape> shape : MyBox3) {
		prims.push_back(std::make_shared<GeometricPrimitive>(shape, myMaterial3, nullptr,
			MediumInterface()));
	}

	std::shared_ptr<Primitive> bvh = std::make_shared<BVHAccel>(prims);

	Transform ltw = Translate(Tangent3f(0, 0, 2.5));
	std::shared_ptr<Light> myLight = std::make_shared<PointLight>(ltw, nullptr, Spectrum(10, 10, 10));
	std::vector<std::shared_ptr<Light>> myLights = { myLight };

	Scene myScene(bvh, myLights);

	std::shared_ptr<Sampler> mySampler = std::make_shared<SobolSampler>(SPP, sampleBounds);

	PathIntegrator myIntegrator(5, myCam, mySampler, Bounds2i({ 0, 0 }, uhd));
	myIntegrator.Render(myScene);

	return 0;
}