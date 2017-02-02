#pragma once

#include <type_traits>
#include <stdexcept>

namespace hypergeom {

// Compile-time Digamma(z) for positive integer argument z
constexpr double psiinteger(int z) {
    return
        (z == 1) ? -0.57721566490153286061 :
        (z > 1) ? psiinteger(z - 1) + 1.0 / (z - 1.0):
        throw std::invalid_argument("Not defined for nonpositive integers");
}

// Compile-time Digamma(z+1/2) for integer argument z
constexpr double psihalfint(int z) {
    return
        (z == 0) ? -1.9635100260214234794410 :
        (z > 0) ? psihalfint(z - 1) + 1.0 / (z - 0.5):
        psihalfint(z + 1) - 1.0 / (z + 0.5);
}

// Compile-time Digamma(z) for positive integer or half-integer argument z2/2
constexpr double psihalf(int z2) {
    return
        (z2 % 2 == 0) ? psiinteger(z2 / 2) :
        psihalfint((z2 - 1) / 2);
}

// Force evaluation of psihalf
template<int twicez>
double psi() {
    constexpr int z2 = std::integral_constant<int, twicez>::value;
    constexpr double ret = psihalf(z2);
    return ret;
}

}
