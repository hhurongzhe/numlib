#pragma once
#ifndef GAUSS_LEGENDRE_HPP
#define GAUSS_LEGENDRE_HPP

#include <iostream>
#include <vector>

namespace gauss_legendre
{
    // Legendre polynomials P(n,x)
    double poly(const int degree, const double x)
    {
        if (degree == 0)
        {
            return 1.0;
        }

        double poly_im1 = 0.0;

        double poly_i = 1.0;

        double poly_ip1 = x;

        for (int i = 1; i < degree; i++)
        {
            const double ip1 = i + 1;

            const double two_ip1 = 2 * i + 1;

            const double one_over_ip1 = 1.0 / ip1;

            poly_im1 = poly_i;

            poly_i = poly_ip1;

            poly_ip1 = (two_ip1 * x * poly_i - i * poly_im1) * one_over_ip1;
        }

        return poly_ip1;
    }

    // Derivative of the Legendre polynomial P'(n,x)
    double poly_der(const int degree, const double x)
    {
        if (degree == 0)
        {
            return 0.0;
        }

        double poly_im1 = 0.0;

        double poly_i = 1.0;

        double poly_ip1 = x;

        double poly_der_i = 0.0;

        double poly_der_ip1 = 1.0;

        for (int i = 1; i < degree; i++)
        {
            const double ip1 = i + 1;

            const double two_ip1 = 2 * i + 1;

            const double one_over_ip1 = 1.0 / ip1;

            poly_im1 = poly_i;

            poly_i = poly_ip1;

            poly_der_i = poly_der_ip1;

            poly_ip1 = (two_ip1 * x * poly_i - i * poly_im1) * one_over_ip1;

            poly_der_ip1 = x * poly_der_i + ip1 * poly_i;
        }

        return poly_der_ip1;
    }

    // associated Legendre polynomials: P_l^1(x)
    double asso_poly1(const int degree, const double x)
    {
        double temp;
        temp = -(sqrt(1.0 - x * x)) * gauss_legendre::poly_der(degree, x);
        return temp;
    }

    // associated Legendre polynomials: P_l^2(x)
    double asso_poly2(const int degree, const double x)
    {
        double temp;
        temp = -degree * gauss_legendre::poly_der(degree + 1, x) + (degree + 2) * x * gauss_legendre::poly_der(degree, x);
        return temp;
    }

} // namespace gauss_legendre

#endif // GAUSS_LEGENDRE_HPP