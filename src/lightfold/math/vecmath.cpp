#include <math/transform.h>
#include <math/vecmath.h>

namespace lightfold {

    // DirectionCone Function Definitions
    DirectionCone Union(const DirectionCone& a, const DirectionCone& b) {
        // Handle the cases where one or both cones are empty
        if (a.IsEmpty())
            return b;
        if (b.IsEmpty())
            return a;

        // Handle the cases where one cone is inside the other
        float theta_a = SafeACos(a.cosTheta), theta_b = SafeACos(b.cosTheta);
        float theta_d = AngleBetween(a.w, b.w);
        if (std::min(theta_d + theta_b, Pi) <= theta_a)
            return a;
        if (std::min(theta_d + theta_a, Pi) <= theta_b)
            return b;

        // Compute the spread angle of the merged cone, $\theta_o$
        float theta_o = (theta_a + theta_d + theta_b) / 2;
        if (theta_o >= Pi)
            return DirectionCone::EntireSphere();

        // Find the merged cone's axis and return cone union
        float theta_r = theta_o - theta_a;
        Tangent3f wr = Cross(a.w, b.w);
        if (LengthSquared(wr) == 0)
            return DirectionCone::EntireSphere();
        Tangent3f w = Rotate(theta_r, wr)(a.w);
        return DirectionCone(w, std::cos(theta_o));
    }

} // namespace lightfold
