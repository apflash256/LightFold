#include <cstdlib>
#include <cmath>

namespace lightfold{

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

} // namespace lightfold