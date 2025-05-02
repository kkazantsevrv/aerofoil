#include <vector>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#pragma once

using mpint = int;

void print(std::vector<double> vec);

void linear_solve (std::vector<double> a,
       std::vector<double> b,
       std::vector<double> c, 
       std::vector<double> d,
       std::vector<double>& x,
       double c1, double c2);

class CubicSpline{
public:
    CubicSpline(std::vector<double> x, std::vector<double> y): _x(x), _y(y){};
    CubicSpline() = default;
    double div(size_t i, size_t j, size_t k);
    void assemble();
    std::vector<double> sinta();
    double sintall();
    double spl(const double& point);
    
    ~CubicSpline() = default;
private:
    std::vector<double> _x;
    std::vector<double> _h;
    std::vector<double> _y;

    std::vector<double> _M;
    std::vector<double> _mu;
    std::vector<double> _lam;
    std::vector<double> _d;

};

