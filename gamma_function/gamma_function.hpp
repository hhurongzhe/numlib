#pragma once
#ifndef GAMMA_FUNCTION_HPP
#define GAMMA_FUNCTION_HPP

#include <cmath>
#include <complex>
#include <iostream>

namespace gamma_function
{
    constexpr double Pi = 3.14159265358979323846;
    //---------------------------------------------------------------
    std::complex<double> log_gamma(const std::complex<double> &z)
    {
        const double x = real(z);
        const double y = imag(z);

        if (x >= 0.5)
        {
            const double log_sqrt_2Pi = 0.91893853320467274177;

            const double g = 4.7421875;

            const std::complex<double> z_m_0p5 = z - 0.5;

            const std::complex<double> z_pg_m0p5 = z_m_0p5 + g;

            const std::complex<double> zm1 = z - 1.0;

            const double c[15] = {0.99999999999999709182, 57.156235665862923517, -59.597960355475491248,
                                  14.136097974741747174, -0.49191381609762019978, 0.33994649984811888699e-4,
                                  0.46523628927048575665e-4, -0.98374475304879564677e-4, 0.15808870322491248884e-3,
                                  -0.21026444172410488319e-3, 0.21743961811521264320e-3, -0.16431810653676389022e-3,
                                  0.84418223983852743293e-4, -0.26190838401581408670e-4, 0.36899182659531622704e-5};

            std::complex<double> sum = c[0];

            for (int i = 1; i < 15; i++)
                sum += c[i] / (zm1 + double(i));

            const std::complex<double> log_Gamma_z = log_sqrt_2Pi + log(sum) + z_m_0p5 * log(z_pg_m0p5) - z_pg_m0p5;

            return log_Gamma_z;
        }
        else if (y >= 0.0)
        {
            const int n = (x < rint(x)) ? (static_cast<int>(rint(x)) - 1) : (static_cast<int>(rint(x)));

            const double log_Pi = 1.1447298858494002;

            const std::complex<double> log_const(-0.693147180559945309417, Pi / 2.0);

            const std::complex<double> i_Pi(0.0, Pi);

            const std::complex<double> eps = z - double(n);

            const std::complex<double> log_sin_Pi_z =
                (y > 110) ? (-i_Pi * z + log_const) : (log(sin(Pi * eps)) - i_Pi * double(n));

            const std::complex<double> log_Gamma_z = log_Pi - log_sin_Pi_z - gamma_function::log_gamma(1.0 - z);

            return log_Gamma_z;
        }
        else
        {
            return conj(gamma_function::log_gamma(conj(z)));
        }
    }
    //---------------------------------------------------------------

    // gamma function.
    std::complex<double> gamma(const std::complex<double> &z) { return exp(gamma_function::log_gamma(z)); }
} // namespace gamma_function

#endif // GAMMA_FUNCTION_HPP