#include "./numlib_define/numlib_define.hpp"

void benchmark()
{
    for (int ix = 0; ix <= 10; ix = ix + 1)
    {
        for (int iy = 0; iy <= 10; iy = iy + 1)
        {
            std::complex<double> z{ix, iy};
            auto value1 = gamma_function::gamma(z);
            std::cout << value1 << std::endl;
        }
    }
}

int main()
{
    benchmark();
    return 0;
}