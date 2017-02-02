#pragma once

#include "gamma.h"
#include "psi.h"

#include <stdexcept>
#include <cmath>

#ifndef DOTRACE
# define DOTRACE 0
#endif

#if DOTRACE
# include <iostream>
# define TRACE do { std::cout << __PRETTY_FUNCTION__ << std::endl; } while (false)
#else
# define TRACE do { } while (false)
#endif

namespace hypergeom {

// Evaluate 2F1(a2/a, b2/2, c2/2, z) using truncated series
template<int a2, int b2, int c2>
double hyp2f1ser(double z) {
    static_assert(c2 > 0 || (c2 % 2 != 0), "c cannot be nonpositive integer");
    static_assert(a2 > 0 || (a2 % 2 != 0), "finite case (check a)");
    static_assert(b2 > 0 || (b2 % 2 != 0), "finite case (check b)");
    TRACE;
    const double a = a2 / 2.;
    const double b = b2 / 2.;
    const double c = c2 / 2.;

    double s = 1;
    double olds = 0;
    double A = a * b / c * z;

    for (int m = 1; olds != s; m++) {
        olds = s;
        s += A;
        A *= (a+m) * (b+m) / (c+m) / (1+m) * z;
    }

    return s;
}

template<int a2, int b2, int n2>
double hyp2f1ser2tail(double w, double wn) {
    TRACE;
    const double a = 0.5 * a2;
    const double b = 0.5 * b2;
    constexpr int n = n2 / 2;
    constexpr double sign = ((n % 2) == 0) ? -1.0 : 1.0;

    double B = sign * wn * gamma<a2 + n2>() * gamma<b2 + n2>()
                / gamma<a2>() / gamma<b2>() / gamma<n2+2>();
    double C = std::log(w) + psi<a2+n2>() + psi<b2+n2>()
                - psi<2>() - psi<2+n2>();

    double s = B * C;
    double olds = 0;
    for (int m = n + 1; s != olds; m++) {
        B *= (a+m-1) * (b+m-1) / m / (m - n) * w;
        C += 1.0 / (a+m-1) + 1.0 / (b+m-1) - 1.0 / (m - n) - 1.0 / m;
        olds = s;
        s += B * C;
    }

    return s;
}

// Evaluate 2F1(a2/2, b2/2; (n2+a2+b2)/2; 1-w) using truncated series
template<int a2, int b2, int n2>
struct hyp2f1ser2;

template<int a2, int b2>
struct hyp2f1ser2<a2, b2, 0> { static double evaluate(double w) {
    TRACE;

    const double mul = gamma<a2+b2>() / gamma<a2>() / gamma<b2>();

    return mul * hyp2f1ser2tail<a2, b2, 0>(w, 1.0);
} };

template<int a2, int b2, int n2>
struct hyp2f1ser2 { static double evaluate(double w) {
    TRACE;
    static_assert(n2 > 0, "n <= 0 in hyp2f1ser2 specialization for n > 0");
    const double a = 0.5 * a2;
    const double b = 0.5 * b2;
    double wn = w;
    constexpr int n = n2 / 2;

    const double mul = gamma<n2+a2+b2>() / gamma<n2+a2>() / gamma<n2+b2>();
    double s = 0;

    double A = gamma<n2>();
    for (int m = 0; m < n-1; m++) {
        s += A;
        wn *= w;
        A = -A * (a+m) * (b+m) / (1+m) / (n-m-1) * w;
    }
    s += A;

    return mul * (s + hyp2f1ser2tail<a2, b2, n2>(w, wn));
} };

// Use reflection method to compute 2F1(a2/2, b2/2; c2/2; 1-w) for w <= 1/2
template<int a2, int b2, int c2, bool k_is_integer, bool k_geq_0>
struct hyp2f1refl;

template<int a2, int b2, int c2, bool k_geq_0>
struct hyp2f1refl<a2, b2, c2, /* k in Z = */ false, k_geq_0> { static double evaluate(double w) {
    TRACE;
    constexpr int k2 = std::integral_constant<int, c2 - a2 - b2>::value;
    const double m1 = gamma<k2>() / gamma<c2-a2>() / gamma<c2-b2>();
    const double m2 = gamma<-k2>() / gamma<a2>() / gamma<b2>() * std::pow(w, 0.5*k2);

    const double f1 = hyp2f1ser<a2, b2, 2 - k2>(w);
    const double f2 = hyp2f1ser<c2-a2, c2-b2, 2 + k2>(w);

    return gamma<c2>() * (m1 * f1 + m2 * f2);
} };

template<int a2, int b2, int c2>
struct hyp2f1refl<a2, b2, c2, /* k in Z */ true, /* k >= 0*/ true> { static double evaluate(double w) {
    TRACE;
    return hyp2f1ser2<a2, b2, c2-a2-b2>::evaluate(w);
} };

template<int a2, int b2, int c2>
struct hyp2f1refl<a2, b2, c2, /* k in Z */ true, /* k >= 0 */ false> { static double evaluate(double w) {
    TRACE;
    constexpr int k = std::integral_constant<int, (c2 - a2 - b2) / 2>::value;
    return std::pow(w, k) * hyp2f1refl<c2-a2, c2-b2, c2, /* k in Z */ true, /* k >= 0 */true>::evaluate(w);
} };

// Evaluate when z in [0, 1)
template<int a2, int b2, int c2>
double hyp2f1roc(double z) {
    TRACE;
    static_assert(b2 > 0, "b should be positive");
    static_assert(a2 > 0, "a should be positive");
    static_assert(c2 > b2, "c should be greater than b");

    if (z < 0.5)
        return hyp2f1ser<a2, b2, c2>(z);

    return hyp2f1refl<a2, b2, c2,
                (c2-a2-b2) % 2 == 0,
                c2 >= a2 + b2
            >::evaluate(1-z);
}

// Dispatch between finite (polynomial) and infinite (series) case
template<int a2, int b2, int c2, bool finite>
struct hyp2f1decide;

template<int a2, int b2, int c2>
struct hyp2f1decide<a2, b2, c2, /* finite = */ true> { static double evaluate(double z) {
    TRACE;
    const double zc = 1-z;

    constexpr int p = std::integral_constant<int, (a2 - c2) / 2>::value;
    const double a = p;
    const double b = 0.5 * b2;
    const double c = 0.5 * c2;
    const double w = z / zc;

    const double mul = std::pow(zc, -b);
    double s = 1;
    double adm = a;
    double bm = b;
    double cm = c;
    double mf = 1;
    double wm = w;

    for (int m = 1; m < p; m++) {
        s += adm * bm / cm / mf * wm;
        adm *= p-m;
        bm *= b+m;
        cm *= c+m;
        mf *= 1+m;
        wm *= w;
    }
    s += adm * bm / cm / mf * wm;

    return mul * s;
} };

template<int a2, int b2, int c2>
struct hyp2f1decide<a2, b2, c2, /* finite = */ false> { static double evaluate(double z) {
    TRACE;
    if (z < 0) // Perform Pfaff transform
        return std::pow(1-z, -0.5*a2) * hyp2f1roc<a2, c2-b2, c2>(z/(z-1));

    return hyp2f1roc<a2, b2, c2>(z);
} };


template<int a2, int b2, int c2>
double hyp2f1(double z) {
    TRACE;
    static_assert(b2 > 0, "b should be positive");
    static_assert(a2 > 0, "a should be positive");
    static_assert(c2 > b2, "c should be greater than b");
    if (z >= 1)
        throw std::invalid_argument("z should be less than 1");

    return hyp2f1decide<a2, b2, c2,
                ((c2 - a2) <= 0) && (((c2 - a2) % 2) == 0)
            >::evaluate(z);
}

}
