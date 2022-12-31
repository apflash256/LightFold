#include <cstdlib>
#include <cmath>

namespace lightfold{

constexpr float Pi = 3.14159265358979323846f;
constexpr float InvPi = 0.31830988618379067154f;
constexpr float Inv2Pi = 0.15915494309189533577f;
constexpr float Inv4Pi = 0.07957747154594766788f;
constexpr float PiOver2 = 1.57079632679489661923f;
constexpr float PiOver4 = 0.78539816339744830961f;
constexpr float Sqrt2 = 1.41421356237309504880f;

static constexpr float Infinity = std::numeric_limits<float>::infinity();
static constexpr float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5f;

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto DifferenceOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto differenceOfProducts = std::fma(a, b, -cd);
    auto error = std::fma(-c, d, cd);
    return differenceOfProducts + error;
}

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto SumOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto sumOfProducts = std::fma(a, b, cd);
    auto error = std::fma(c, d, -cd);
    return sumOfProducts + error;
}

template <typename T, typename U, typename V>
inline constexpr T Clamp(T val, U low, V high) {
    if (val < low)
        return T(low);
    else if (val > high)
        return T(high);
    else
        return val;
}
// http://www.plunk.org/~hatch/rightway.html
inline float SinXOverX(float x) {
    if (1 + x * x == 1)
        return 1;
    return std::sin(x) / x;
}
inline float SafeASin(float x) {
    return std::asin(Clamp(x, -1, 1));
}
inline float SafeACos(float x) {
    return std::acos(Clamp(x, -1, 1));
}
inline double SafeASin(double x) {
    return std::asin(Clamp(x, -1, 1));
}
inline double SafeACos(double x) {
    return std::acos(Clamp(x, -1, 1));
}

inline float SafeSqrt(float x) {
    return std::sqrt(std::max(0.f, x));
}
inline double SafeSqrt(double x) {
    return std::sqrt(std::max(0., x));
}

inline constexpr float gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}


} // namespace lightfold