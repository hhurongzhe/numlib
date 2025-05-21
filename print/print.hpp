#pragma once
#ifndef PRINT_HPP
#define PRINT_HPP

#include <iostream>

namespace print
{
    void print_vector(std::vector<double> list)
    {
        int num = list.size();
        for (auto x : list)
        {
            std::cout << x << std::endl;
        }
    }

} // namespace print

#endif // PRINT_HPP