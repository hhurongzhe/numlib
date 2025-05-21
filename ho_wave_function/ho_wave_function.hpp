#pragma once
#ifndef HO_WF_HPP
#define HO_WF_HPP

#include "./../gamma_function/gamma_function.hpp"
#include "./../basic_math/basic_math.hpp"
#include "./../gauss_laguerre/gauss_laguerre.hpp"
#include <iostream>

namespace ho_wave_function
{
    using basic_math::iphase;
    using std::exp;
    using std::pow;
    using std::sqrt;

    double ho_wave_function_momentumspace(const int n, const int l, const double length, const double q)
    {
        double temp;
        const double b = length;
        const double bq = length * q;
        const double bqbq = bq * bq;
        constexpr double sqrt_pi = 1.772453850905516;
        const double num = basic_math::factorial(n) * basic_math::factorial(n + l);
        const double den = sqrt_pi * basic_math::factorial(2 * l + 2 * n + 1);
        const double part_1 = sqrt(pow(b, 3)) * basic_math::pow_two(n + l + 1);
        const double part_2 = sqrt(num / den);
        const double part_3 = pow(bq, l);
        const double part_4 = exp(-0.5 * bqbq);
        const double part_5 = gauss_laguerre::poly(n, l + 0.5, bqbq);
        temp = part_1 * part_2 * part_3 * part_4 * part_5;
        return temp;
    }

    double ho_wave_function_coordinatespace(const int n, const int l, const double alpha, const double r)
    {
        double temp;
        //* a = 1/b = sqrt(mu*Omega).
        const double a = alpha;
        const double ar = alpha * r;
        const double arar = ar * ar;
        constexpr double sqrt_pi = 1.772453850905516;
        const double num = basic_math::factorial(n) * basic_math::factorial(n + l);
        const double den = sqrt_pi * basic_math::factorial(2 * l + 2 * n + 1);
        const double part_1 = sqrt(pow(a, 3)) * basic_math::pow_two(n + l + 1);
        const double part_2 = sqrt(num / den);
        const double part_3 = pow(ar, l);
        const double part_4 = exp(-0.5 * arar);
        const double part_5 = gauss_laguerre::poly(n, l + 0.5, arar);
        temp = part_1 * part_2 * part_3 * part_4 * part_5;
        return temp;
    }

} // namespace HO_wave_function

#endif // HO_WF_HPP