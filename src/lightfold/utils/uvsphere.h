#pragma once

#include <core/shape.h>
#include <core/light.h>

namespace lightfold {

	std::vector<std::shared_ptr<Primitive>> UVSphere(const Transform* ObjectToWorld,
		const Transform* WorldToObject, int lats, int longs,
		const std::shared_ptr<Material>& material,
		const std::shared_ptr<AreaLight>& areaLight,
		const MediumInterface& mediumInterface);

} // namespace lightfold