#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>
#include <tuple>
#include "spline.hpp"
#include "utils.hpp"

constexpr double pi = acos(-1.0);

std::vector<double> set_grid(size_t _n_points, double l, double r){
    std::vector<double> _s(_n_points, 0.0);
    for(size_t i=0; i<_n_points; i++){
        _s[i] = l + i*(r-l)/(_n_points-1.0);
    }
    return _s;
}
class Aerofoil{
public:
    Aerofoil(size_t n, double l=0.0, double r=1.0);
    void set_v();
    void compute_fi();
    void compute_params();
    void compute_figam();
private:
    size_t _n_points;
    double _sa;
    std::vector<double> _s;
    std::vector<double> _gam;
    std::vector<double> _v;
    std::vector<double> _fi;
    std::vector<double> _figam;

    CubicSpline _vs_spl;
    
    double _Gamma;
    double _beta;
    double _gamma_a;
    double _u0;
    double _c1;
};

Aerofoil::Aerofoil(size_t n, double l, double r){
    _n_points = n;
    _s = std::vector<double>(_n_points, 0.0);
    _v = std::vector<double>(_n_points, 0.0);
    _fi = std::vector<double>(_n_points, 0.0);
    _figam = std::vector<double>(_n_points, 0.0);
    _s = set_grid(_n_points, l, r);
    _gam = set_grid(_n_points, 0.0, 2.0*pi);
    for(size_t i=0; i<_n_points; i++){
        _v[i] = _s[i]*_s[i];
    }
    for(size_t i=0; i<_n_points-1; i++){
        if(_v[i]*_v[i+1] < 0){
            _sa = _s[i] - _v[i]/(_v[i+1]-_v[i])*(_s[i+1]-_s[i]);
            break;
        }
    }
    _vs_spl = CubicSpline(_s, _v);
    _vs_spl.assemble();
}

void Aerofoil::compute_fi(){
    std::vector<double> s = set_grid(_n_points, 0.0, _sa);
    std::vector<double> v(_n_points, 0.0);
    for(size_t i=0; i<_n_points; i++){
        v[i] = _vs_spl.spl(s[i]);
    }
    auto spl = CubicSpline(s, v);
    double fi0 = spl.sintall();
    spl = CubicSpline(_s, _v);
    _fi = spl.sinta();
    for(size_t i=0; i<_fi.size(); i++){
        _fi[i]-= fi0;
    }
}

void Aerofoil::compute_params(){
    double _Gamma = _fi[_n_points-1] - _fi[0];
    double a = pi/_Gamma*_fi[_n_points-1]-pi/2.0;
    if(_Gamma > 0){
        double zero = pi/4.0;
        //_beta = newton(zero, 1e-12, 10, a);
        _beta = newton(
                std::function<double(double, double)>(
                    [](double beta, double param) { return func(beta, param);}
                ),
                std::function<double(double, double)>(
                    [](double beta, double param) { return dfunc(beta, param);}
                ),
                zero, 1e-12, 10, a);
    }else if(_Gamma < 0){
        double zero = -pi/4.0;
        //_beta = newton(zero, 1e-12, 10, a);
        _beta = newton(
                std::function<double(double, double)>(
                    [](double beta, double param) { return func(beta, param);}
                ),
                std::function<double(double, double)>(
                    [](double beta, double param) { return dfunc(beta, param);}
                ),
                zero, 1e-12, 10, a);
    }else if(_Gamma - 0.0 < 1e-12){
        _beta = 0.0;
        _gamma_a = 0.0;
        _c1 = _fi[_n_points-1]/2.0;
        _u0 = _fi[_n_points-1]/4.0;
        return;
    }
    _gamma_a = pi + 2.0*_beta;
    _u0 = _Gamma/(4.0*pi*sin(_beta));
    _c1 = _fi[_n_points-1] - 2.0*_u0*cos(_beta);
}

void Aerofoil::compute_figam(){
    for(size_t i=0; i<_n_points; i++){
        double zero;
        if(_s[i] < _sa){
            zero = _sa/2.0;
        }else{
            zero = _s[_n_points-1] - _sa/2.0;
        }
        _figam[i] =  newton<double, double, double, double, double, double>(
                std::function<double(double, double, double, double, double, double)>(
                    [](double gam, double beta, double u0, double Gamma, double c1, double fi) 
                    { return func1(gam, beta, u0, Gamma, c1, fi);}
                ),
                std::function<double(double, double, double, double, double, double)>(
                    [](double gam, double beta, double u0, double Gamma, double c1, double fi) 
                    { return dfunc1(gam, beta, u0, Gamma, c1, fi);}
                ),
                zero, 1e-12, 10, _gam[i], _beta, _u0, _Gamma, _c1, _fi[i]);
    }
}


int main(){
   std::cout << " aaaaa " << std::endl;
}
