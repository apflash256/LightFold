#pragma once

#include <cstdint>

namespace lightfold {

	static constexpr int NumSobolDimensions = 1024;
	static constexpr int SobolMatrixSize = 52;
	extern const uint32_t SobolMatrices32[NumSobolDimensions * SobolMatrixSize];
	extern const uint64_t SobolMatrices64[NumSobolDimensions * SobolMatrixSize];
	extern const uint64_t VdCSobolMatrices[][SobolMatrixSize];
	extern const uint64_t VdCSobolMatricesInv[][SobolMatrixSize];

} // namespace lightfold