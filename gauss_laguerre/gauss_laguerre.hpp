#pragma once
#ifndef GAUSS_LAGUERRE_HPP
#define GAUSS_LAGUERRE_HPP

#include <iostream>
#include <vector>

namespace gauss_laguerre
{

    // Calculation of the Laguerre polynomial L(n, alpha, x)
    // ---------------------------------------------- -
    // The following recurrence relation is used :
    // (i+1).L(i+1, alpha, x) = (-x + alpha + 1 + 2i).L(i, alpha, x) - (i+a).L(i-1, alpha, x) for i > 0.
    // L(0, alpha, x) = 1. This value is directly returned if the degree of the polynomial is 0.
    // L(1, alpha, x) = alpha + 1 - x. This value is directly returned if the degree of the polynomial is 1.
    //
    //
    // Variables:
    // ----------
    // degree : degree of the Laguerre polynomial
    // i, poly_im1, poly_i, poly_ip1 : Laguerre polynomials of degree i-1, i, and i+1, with n > 1.
    // i goes from 1 to n - 1.
    // alpha : parameter of the Laguerre polynomial
    // minus_x_plus_alpha_plus_one : -x + alpha + 1
    // x : variable of the Laguerre polynomial.

    double poly(const int degree, const double alpha, const double x)
    {
        if (degree == 0)
            return 1.0;

        const double minus_x_plus_alpha_plus_one = -x + alpha + 1.0;

        double poly_im1 = 0.0;

        double poly_i = 1.0;

        double poly_ip1 = minus_x_plus_alpha_plus_one;

        for (int i = 1; i < degree; i++)
        {
            poly_im1 = poly_i;

            poly_i = poly_ip1;

            poly_ip1 = ((minus_x_plus_alpha_plus_one + 2 * i) * poly_i - (i + alpha) * poly_im1) / (i + 1.0);
        }

        return poly_ip1;
    }

    // Calculation of the Laguerre polynomial derivative dL/dx(n, alpha, x)
    // --------------------------------------------------------------
    // The following recurrence relations are used :
    // (i+1).L'(i+1, alpha, x) = (-x + alpha + 1 + 2i).L'(i, alpha, x) - (i+a).L'(i-1, alpha, x) - L(i, alpha, x) for i > 0.
    // (i+1).L(i+1, alpha, x) = (-x + alpha + 1 + 2i).L(i, alpha, x) - (i+a).L(i-1, alpha, x) for i > 0.
    // L(0, alpha, x) = 1
    // L(1, alpha, x) = 1 + alpha - x
    // L'(0, alpha, x) =  0. This value is directly returned if the degree of the polynomial is 0.
    // L'(1, alpha, x) = -1. This value is directly returned if the degree of the polynomial is 1.
    //
    // I do not use the formula L'(i+1 , alpha , x) = ((i+1).L(i+1 , alpha , x) - (i + 1 + alpha).L(i , alpha , x))/x as it
    // is unstable for x ~ 0.
    //
    // Variables:
    // ----------
    // degree : degree of the Laguerre polynomial, not of the derivative of the Laguerre polynomial.
    // i, poly_im1, poly_i, poly_ip1, poly_der_im1, poly_der_i, poly_der_ip1 : Laguerre polynomials and derivatives of
    // degree i-1, i, and i+1, with n > 1. i goes from 1 to n - 1. one_over_ip1 : 1/(i+1) minus_x_plus_alpha_plus_one : -x +
    // alpha + 1 minus_x_plus_alpha_plus_one_plus_2i : -x + alpha + 1 + 2i x : variable of the Laguerre polynomial
    // derivative. alpha : parameter of the Laguerre polynomial

    double poly_der(const int degree, const double alpha, const double x)
    {
        if (degree == 0)
            return 0.0;

        const double minus_x_plus_alpha_plus_one = -x + alpha + 1.0;

        double poly_im1 = 0.0;

        double poly_i = 1.0;

        double poly_ip1 = minus_x_plus_alpha_plus_one;

        double poly_der_im1 = 0.0;

        double poly_der_i = 0.0;

        double poly_der_ip1 = -1.0;

        for (int i = 1; i < degree; i++)
        {
            const double one_over_ip1 = 1.0 / (i + 1.0);

            const double i_plus_alpha = i + alpha;

            const double minus_x_plus_alpha_plus_one_plus_2i = minus_x_plus_alpha_plus_one + 2.0 * i;

            poly_im1 = poly_i;

            poly_der_im1 = poly_der_i;

            poly_i = poly_ip1;

            poly_der_i = poly_der_ip1;

            poly_ip1 = (minus_x_plus_alpha_plus_one_plus_2i * poly_i - i_plus_alpha * poly_im1) * one_over_ip1;

            poly_der_ip1 =
                (minus_x_plus_alpha_plus_one_plus_2i * poly_der_i - i_plus_alpha * poly_der_im1 - poly_i) * one_over_ip1;
        }

        return poly_der_ip1;
    }

} // namespace gauss_laguerre

#endif // GAUSS_LAGUERRE_HPP