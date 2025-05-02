#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>
#include <tuple>
#include "spline.hpp"
#include "utils.hpp"
#include "plot.hpp"
#include <fstream>
#include <string>
constexpr double pi = acos(-1.0);

void print_f(std::string filename, std::vector<double> x, std::vector<double> y){
    std::ofstream out;
    out.open(filename);
    for(size_t i=0; i<x.size(); i++){
        out << x[i] << "\t" << y[i] << "\n";
    }
    out.close();
}

std::vector<double> set_v(const std::vector<double>& s){
    auto v = std::vector<double>(s.size(), 0.0);
    for(size_t i=0; i<s.size(); i++){
        v[i] = -sin(5.0*s[i]) + sin(5.0)/2.0;
    }
    return v;
}

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
    void step();
private:
    size_t _n_points;
    double _sa;
    std::vector<double> _s;
    std::vector<double> _gam;
    std::vector<double> _v;
    std::vector<double> _fi;
    std::vector<double> _gams;
    std::vector<double> _vgam;

    CubicSpline _vs_spl;
    CubicSpline _vgam_spl;
    
    double _Gamma;
    double _beta;
    double _gamma_a;
    double _u0;
    double _c1;
    
    void compute_fi();
    void compute_params();
    void compute_gams();
    void compute_vgam();
};

Aerofoil::Aerofoil(size_t n, double l, double r){
    _n_points = n;
    _s = std::vector<double>(_n_points, 0.0);
    _v = std::vector<double>(_n_points, 0.0);
    _fi = std::vector<double>(_n_points, 0.0);
    _gams = std::vector<double>(_n_points, 0.0);
    _vgam = std::vector<double>(_n_points, 0.0);
    _s = set_grid(_n_points, l, r);
    _gam = set_grid(_n_points, 0.0, 2.0*pi);
    _v = set_v(_s);
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
    spl.assemble();
    double fi0 = spl.sintall();
    spl = CubicSpline(_s, _v);
    spl.assemble();
    _fi = spl.sinta();
    for(size_t i=0; i<_fi.size(); i++){
        _fi[i]-= fi0;
    }
}

void Aerofoil::compute_params(){
    _Gamma = _fi[_n_points-1] - _fi[0];
    double a = pi/_Gamma*_fi[_n_points-1]-pi/2.0;
    struct p params;
    params.param = a;
    if(_Gamma > 0){
        double zero = pi/4.0;
        _beta = newton(
                std::function<double(double, struct p)>(
                    [](double beta, struct p p) { return func(beta, p);}
                ),
                std::function<double(double, struct p)>(
                    [](double beta, struct p p) { return dfunc(beta, p);}
                ),
                zero, 1e-12, 10, params
                );
    }else if(_Gamma < 0){
        double zero = -pi/4.0;
        _beta = newton(
                std::function<double(double, struct p)>(
                    [](double beta, struct p p) { return func(beta, p);}
                ),
                std::function<double(double, struct p)>(
                    [](double beta, struct p p) { return dfunc(beta, p);}
                ),
                zero, 1e-12, 10, params
                );
    }else if(_Gamma - 0.0 < 1e-12){
        _beta = 0.0;
        _gamma_a = pi;
        _c1 = _fi[_n_points-1]/2.0;
        _u0 = _fi[_n_points-1]/4.0;
        return;
    }
    _gamma_a = pi + 2.0*_beta;
    _u0 = _Gamma/(4.0*pi*sin(_beta));
    _c1 = _fi[_n_points-1] - 2.0*_u0*cos(_beta);
}

void Aerofoil::compute_gams(){
    struct p params;
    params.beta = _beta;
    params.u0 = _u0;
    params.c1 = _c1;
    params.Gamma = _Gamma;
    for(size_t i=0; i<_n_points; i++){
        double zero;
        if(_gam[i] < _gamma_a){
            zero = _gamma_a/2.0;
        }else{
            zero = (2.0*pi + _gamma_a)/2.0;
        }
        params.fi = _fi[i];
        if(_s[i] < _sa) zero = pi + _gamma_a/2.0;
        else zero = _gamma_a/2.0;
        _gams[i] =  newton(
                std::function<double(double, struct p)>(
                    [](double gam, struct p p) 
                    { return func1(gam, p);}
                    ),
                std::function<double(double, struct p)>(
                    [](double gam, struct p p) 
                    { return dfunc1(gam, p);}
                    ),
                zero, 1e-12, 20, params);
    }
    _gams[0] = 2.0*pi;
    _gams[_n_points-1] = 0.0;
    //print_f("dat.dat", _gam, _gams);

}

void Aerofoil::compute_vgam(){
    auto gams_ = _gams;
    auto s_ = _s;
    for(size_t i=0; i< _n_points/2; i++){
        std::swap(gams_[i], gams_[_n_points-i-1]);
        std::swap(s_[i], s_[_n_points-i-1]);
    }
    auto spl = CubicSpline(gams_, s_);
    spl.assemble();
    std::vector<double> sgam(_n_points, 0.0);
    for(size_t i=0; i<_n_points; i++){
        sgam[i] = spl.spl(_gam[i]);
        _vgam[i] = _vs_spl.spl(sgam[i]);
    }
    _vgam_spl = CubicSpline(_gam, _vgam);
    _vgam_spl.assemble();
    //print_f("dat.dat", _gam, _vgam);
}

void Aerofoil::step(){
    compute_fi();
    compute_params();
    compute_gams();
    compute_vgam();
}


int main(){
   std::cout << " aaaaa " << std::endl;
   Aerofoil aero_test(200, 0, 1);
   aero_test.step();
}
