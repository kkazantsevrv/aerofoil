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
#include <complex>

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
    Aerofoil(double vinf, size_t n, double l=0.0, double r=1.0);
    void step();
private:
    size_t _n_points;
    double _vinf;
    double _sa;
    std::vector<double> _s;
    std::vector<double> _gam;
    std::vector<double> _v;
    std::vector<double> _fi;
    std::vector<double> _gams;
    std::vector<double> _vgam;
    std::vector<double> _S_wave;
    std::vector<double> _S1_wave;
    std::vector<double> _theta;
    std::vector<double> _x;
    std::vector<double> _y;

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
    void compute_S();
    void compute_S1();
    void compute_theta();
    void compute_xy();
};

Aerofoil::Aerofoil(double vinf, size_t n, double l, double r): _vinf(vinf){
    _n_points = n;

    _s = std::vector<double>(_n_points, 0.0);
    _v = std::vector<double>(_n_points, 0.0);
    _fi = std::vector<double>(_n_points, 0.0);
    _gams = std::vector<double>(_n_points, 0.0);
    _vgam = std::vector<double>(_n_points, 0.0);
    _S_wave = std::vector<double>(_n_points, 0.0);
    _S1_wave = std::vector<double>(_n_points, 0.0);
    _theta = std::vector<double>(_n_points, 0.0);
    _x = std::vector<double>(_n_points, 0.0);
    _y = std::vector<double>(_n_points, 0.0);

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

void Aerofoil::compute_S(){
    for(size_t i=0; i<_n_points; i++){
        _S_wave[i] = log(std::abs(_vgam[i])) - log(std::abs(2.0*sin((_gam[i] - _gamma_a)/2.0)));
    }
    for(size_t i=1; i<_n_points-1; i++){
        if(std::abs(_gam[i]-_gamma_a) < 1e-9){
            _S_wave[i] = (_S_wave[i-1] + _S_wave[i+1])/2.0;
        }
    }
    //print_f("dat.dat", _gam, _S_wave);
}

void Aerofoil::compute_S1(){
    double A1 = 2.0*pi*log(_vinf);
    double A2 = -pi;
    double A3 = 0.0;
    double mu1;
    double mu2;
    double mu3;

    std::vector<double> S(_S_wave);
    auto spl = CubicSpline(_gam, S);
    spl.assemble();
    mu1 = 1.0/2.0/pi*(spl.sintall() - A1);
    for(size_t i=0; i<_n_points; i++){
        S[i] = _S_wave[i]*cos(_gam[i]);
    }
    spl = CubicSpline(_gam, S);
    spl.assemble();
    mu2 = 1.0/pi*(spl.sintall() - A2);
    for(size_t i=0; i<_n_points; i++){
        S[i] = _S_wave[i]*sin(_gam[i]);
    }
    spl = CubicSpline(_gam, S);
    spl.assemble();
    mu3 = 1.0/pi*(spl.sintall() - A3);

    for(size_t i=0; i<_n_points; i++){
        _S1_wave[i] = _S_wave[i] - (mu1 + mu2*cos(_gam[i]) + mu3*sin(_gam[i]));
    }
}

void Aerofoil::compute_theta(){
    _theta = hilbert(_S1_wave, _gam);
    //print_f("dat.dat", _gam, _theta);
}

void Aerofoil::compute_xy(){
    std::vector<std::complex<double>> f1(_n_points);
    for(size_t i=0; i<_n_points; i++){
        const std::complex<double> I(0.0, 1.0);
        f1[i] = I * _u0
                * std::exp(-I * _beta)
                * (std::exp(I * _gam[i]) - 1.0)
                * std::exp(-_S1_wave[i] - I * _theta[i]);
    }
    std::vector<double> ref1(_n_points);
    std::vector<double> imf1(_n_points);
    for(size_t i=0; i<_n_points; i++){
        ref1[i] = f1[i].real();
        imf1[i] = f1[i].imag();
    }
    auto spl = CubicSpline(_gam, ref1);
    spl.assemble();
    _x = spl.sinta();
    spl = CubicSpline(_gam, imf1);
    spl.assemble();
    _y = spl.sinta();
    print_f("dat.dat", _x, _y);
}

void Aerofoil::step(){
    compute_fi();
    compute_params();
    compute_gams();
    compute_vgam();
    compute_S();
    compute_S1();
    compute_theta();
    compute_xy();
}

int main(){
   std::cout << " aaaaa " << std::endl;
   Aerofoil aero_test(1.0, 1000, 0, 1);
   aero_test.step();
}
