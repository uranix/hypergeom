#pragma once

#include <type_traits>
#include <stdexcept>

namespace hypergeom {

// Compile-time Gamma(z) for positive integer argument z
constexpr double gammainteger(int z) {
    return
        (z == 1) ? 1.0 :
        (z > 1) ? (z - 1) * gammainteger(z - 1) :
        throw std::invalid_argument("Not defined for nonpositive integers");
}

// Compile-time Gamma(z+1/2) for integer argument z
constexpr double gammahalfint(int z) {
    return
        (z == 0) ? 1.7724538509055160273 :
        (z > 0) ? (z - 0.5) * gammahalfint(z - 1) :
        (z % 2 == 0) ? 3.1415926535897932385 / gammahalfint(-z) :
        -3.1415926535897932385 / gammahalfint(-z);
}

// Compile-time Gamma(z) for positive integer or half-integer argument z2/2
constexpr double gammahalf(int z2) {
    return
        (z2 % 2 == 0) ? gammainteger(z2 / 2) :
        gammahalfint((z2 - 1) / 2);
}

// Force evaluation of gammahalf
template<int twicez>
double gamma() {
    constexpr int z2 = std::integral_constant<int, twicez>::value;
    constexpr double ret = gammahalf(z2);
    return ret;
}

}
