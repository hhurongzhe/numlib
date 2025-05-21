#pragma once
#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <cmath>
#include <iostream>

namespace basic_math
{
    constexpr unsigned long long factorial_table[] = {1,
                                                      1,
                                                      2,
                                                      6,
                                                      24,
                                                      120,
                                                      720,
                                                      5040,
                                                      40320,
                                                      362880,
                                                      3628800,
                                                      39916800,
                                                      479001600,
                                                      6227020800,
                                                      87178291200,
                                                      1307674368000,
                                                      20922789888000,
                                                      355687428096000,
                                                      6402373705728000,
                                                      121645100408832000,
                                                      2432902008176640000};

    constexpr unsigned long long doublefactorial_table[] = {1,
                                                            1,
                                                            2,
                                                            3,
                                                            8,
                                                            15,
                                                            48,
                                                            105,
                                                            384,
                                                            945,
                                                            3840,
                                                            10395,
                                                            46080,
                                                            135135,
                                                            645120,
                                                            2027025,
                                                            10321920,
                                                            34459425,
                                                            185794560,
                                                            654729075,
                                                            3715891200,
                                                            13749310575,
                                                            81749606400,
                                                            316234143225,
                                                            1961990553600,
                                                            7905853580625,
                                                            51011754393600,
                                                            213458046676875,
                                                            1428329123020800,
                                                            6190283353629375,
                                                            42849873690624000};

    constexpr double precomputed_sqrt[] = {1.0,
                                           1.41421356237309504880,
                                           1.73205080756887729352,
                                           2.0,
                                           2.23606797749978969640,
                                           2.44948974278317809819,
                                           2.64575131106459059050,
                                           2.82842712474619009760,
                                           3.0,
                                           3.16227766016837933199,
                                           3.31662479035539984911,
                                           3.46410161513775458705,
                                           3.60555127546398929311,
                                           3.74165738677394138558,
                                           3.87298334620741688517,
                                           4.0,
                                           4.12310562561766054982,
                                           4.24264068711928514640,
                                           4.35889894354067355223,
                                           4.47213595499957939281,
                                           4.58257569495584000658,
                                           4.69041575982342955456};

    // factorial function n!, quick for n<=20.
    unsigned long long factorial(unsigned int n)
    {
        if (n <= 20)
        {
            return factorial_table[n];
        }
        return n * factorial(n - 1);
    }

    // double-factorial function n!!, quick for n<=30.
    unsigned long long doublefactorial(unsigned int n)
    {
        if (n <= 30)
        {
            return doublefactorial_table[n];
        }
        return n * doublefactorial(n - 2);
    }

    // quick pow(2,n) function
    unsigned long long pow_two(unsigned int n) { return (1ULL << n); }

    // double -> int
    int make_int(const double n) { return static_cast<int>(rint(n)); }

    // (-1)^n
    double iphase(const int n) { return ((n & 0x01) ? (-1.0) : (1.0)); }

    // sign function
    double sign(const double x) { return ((x < 0) ? (-1.0) : (1.0)); }

    // hat(j) = sqrt(2j+1),
    // which is quick for j <= 21/2.
    double hat(const double j)
    {
        const int two_j = static_cast<int>(2 * j);
        if (two_j <= 21)
        {
            return precomputed_sqrt[two_j];
        }
        else
        {
            return sqrt(2 * j + 1);
        }
    }

    // 计算组合数 C(n, m)
    inline int64_t get_combination_number_C(const int n, const int m)
    {
        if (m > n || m < 0)
            return 0;
        if (m == 0 || m == n)
            return 1;
        std::vector<int64_t> dp(m + 1, 0);
        dp[0] = 1;
        for (int i = 1; i <= n; ++i)
        {
            for (int j = std::min(i, m); j > 0; --j)
            {
                dp[j] = dp[j] + dp[j - 1];
            }
        }
        return dp[m];
    }

} // namespace basic_math

#endif // BASIC_MATH_HPP