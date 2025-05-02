#include <vector>
#include <iostream>
#include "spline.hpp"

int main(){
    size_t n = 10;
    std::vector<double> x = {0, 1, 3, 5, 10};
    std::vector<double> y = {-2, 3, 4, 6, 9};
    auto spl = CubicSpline(x, y);
    spl.assemble();
    print(spl.sinta());
}
