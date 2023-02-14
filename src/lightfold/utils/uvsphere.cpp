#include <utils/uvsphere.h>

#include <core/light.h>
#include <core/primitive.h>
#include <shape/triangle.h>

namespace lightfold {

	std::vector<std::shared_ptr<Primitive>> UVSphere(const Transform* ObjectToWorld,
		const Transform* WorldToObject, int lats, int longs,
		const std::shared_ptr<Material>& material,
		const std::shared_ptr<AreaLight>& areaLight,
		const MediumInterface& mediumInterface) {
		Transform otw = Scale(1, 1, 1), wto = Scale(1, 1, 1);
		Point3f* p = new Point3f[(lats + 1) * longs];
		Normal3f* n = new Normal3f[(lats + 1) * longs];
		Point2f* uv = new Point2f[(lats + 1) * longs];
		for (int i = 0; i < lats + 1; i++) {
			for (int j = 0; j < longs; j++) {
				float theta = Pi * i / lats;
				float phi = 2 * Pi * j / longs;
				float x = std::sin(theta) * std::cos(phi);
				float y = std::sin(theta) * std::sin(phi);
				float z = std::cos(theta);
				p[i * longs + j] = { x, y, z };
				n[i * longs + j] = { x, y, z };
				uv[i * longs + j] = { 1 - (float)i / lats, (float)j / longs };
			}
		}
		int* vInds = new int[6 * (lats - 1) * longs];
		for (int i = 1; i < lats; i++) {
			for (int j = 0; j < longs - 1; j++) {
				vInds[6 * ((i - 1) * longs + j)] = (i - 1) * longs + j;
				vInds[6 * ((i - 1) * longs + j) + 1] = i * longs + j;
				vInds[6 * ((i - 1) * longs + j) + 2] = i * longs + j + 1;
				vInds[6 * ((i - 1) * longs + j) + 3] = (i + 1) * longs + j + 1;
				vInds[6 * ((i - 1) * longs + j) + 4] = i * longs + j + 1;
				vInds[6 * ((i - 1) * longs + j) + 5] = i * longs + j;
			}
			vInds[6 * i * longs - 6] = i * longs - 1;
			vInds[6 * i * longs - 5] = (i + 1) * longs - 1;
			vInds[6 * i * longs - 4] = i * longs;
			vInds[6 * i * longs - 3] = (i + 1) * longs;
			vInds[6 * i * longs - 2] = i * longs;
			vInds[6 * i * longs - 1] = (i + 1) * longs - 1;
		}
		std::vector<std::shared_ptr<Shape>> MyTri = CreateTriangleMesh(
			ObjectToWorld, WorldToObject, false, 2 * (lats - 1) * longs, vInds, (lats + 1) * longs, p, nullptr, n, uv);
		std::vector<std::shared_ptr<Primitive>> prims;
		for (std::shared_ptr<Shape> shape : MyTri) {
			prims.push_back(std::make_shared<GeometricPrimitive>(shape, material, areaLight, mediumInterface));
		}
		return prims;
	}
	
} // namespace lightfold