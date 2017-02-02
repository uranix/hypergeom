#define DOTRACE 1
#include "hyp2f1.h"

#include <iostream>
#include <cmath>
#include <fenv.h>

int main() {
    feenableexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
    std::cout.precision(16);
    std::cout << hypergeom::hyp2f1<4,3,10>(-1) << std::endl;

    return 0;
}
