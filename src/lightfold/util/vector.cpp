#include <util/vector.h>

namespace lightfold {

	Vector Transform(Vector v, Transformation t) {
		Vector res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				res[i] += t.matrix[i][j] * v[j];
		return res;
	}

	Normal Transform(Normal n, Transformation t) {
		Normal res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				res[i] += t.invmat[i][j] * n[j];
		return res;
	}

} // namespace lightfold