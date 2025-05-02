#include <vector>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "spline.hpp"

void print(std::vector<double> vec){
    for(size_t i=0; i<vec.size(); i++)
    std::cout << std::setw(13) << vec[i] << " ";
    std::cout << std::endl;
}

void linear_solve (std::vector<double> a,
       std::vector<double> b,
       std::vector<double> c, 
       std::vector<double> d,
       std::vector<double>& x,
       double c1, double c2){

    
    std::vector<double> aa(a.begin()+1, a.end()-1);
    std::vector<double> bb(b.begin()+1, b.end()-1);
    std::vector<double> cc(c.begin()+1, c.end()-1);
    std::vector<double> dd(d.begin()+1, d.end()-1);
    size_t n = dd.size();

    for(size_t i=1; i<n; i++){
        dd[i]-= dd[i-1]/bb[i-1]*aa[i-1];
        bb[i]-= cc[i-1]/bb[i-1]*aa[i-1];
    }

    for(int i=n-2; i>=0; --i){
        dd[i]-= dd[i+1]/bb[i+1]*cc[i];
    }
    for(size_t i=1; i<n+1; i++){
        x[i] = dd[i-1]/bb[i-1];
    }
    x[0] = (d[0] - c[0]*x[1] - c1*x[2])/b[0];
    x[n+1] = (d[n+1] - a[n]*x[n] - c2*x[n-1])/b[n+1];
}
double CubicSpline::div(size_t i, size_t j, size_t k){
    return (_y[k] - _y[j])/(_x[k] - _x[j])/(_x[k] - _x[i]) - (_y[j] - _y[i])/(_x[j] - _x[i])/(_x[k] - _x[i]);
}

void CubicSpline::assemble(){
    _mu.resize(_x.size() - 1, 0);
    _lam.resize(_x.size() - 1, 0);
    _h.resize(_x.size() - 1, 0);
    _d.resize(_x.size(), 0);
    _M.resize(_x.size(), 0);

    std::vector<double> b(_x.size(), (double)2.0);
    
    for(size_t i=0; i<_x.size()-1; i++){
        _h[i] = (_x[i+1] - _x[i]);
    }

    for(size_t i=1; i<_x.size()-1; i++){
        _mu[i-1] = _h[i-1]/(_h[i-1] + _h[i]);
        _lam[i] = (double)1.0 - _mu[i-1];
        _d[i] = double(6.0)*div(i-1, i, i+1);
    }
    double mu1 = _mu[0];
    double lam1 = _lam[1];
    double mus1 = _mu[_x.size() - 3];
    double lams1 = _lam[_x.size() - 2];
    _mu[_x.size() - 2] = (double)(-1.0);
    _lam[0] = (double)(-1.0);
    _d[0] = (double)0.0;
    _d[_x.size() - 1] = (double)0.0;
    _d[1]*= lam1;
    _d[_x.size() - 2]*= mus1;
    _mu[_x.size() - 3] = mus1 - lams1;
    _mu[0] = (double)0.0;
    _lam[_x.size() - 2] = (double)0.0;
    _lam[1] = lam1 - mu1;
    b[0] = lam1;
    b[1] = (double)1.0 + lam1;
    b[_x.size() - 1] = mus1;
    b[_x.size() - 2] = (double)1.0 + mus1;
    linear_solve(_mu, b, _lam, _d, _M, mu1, lams1);
}

std::vector<double> CubicSpline::sinta(){
    std::vector<double> ls(_x.size(), 0);
    for(size_t i=0; i< _x.size()-1; i++){
        ls[i+1] = -_h[i]*(-2.0*_h[i]*_h[i]*(_M[i+1]+_M[i]) + (_M[i+1]+_M[i])*(_x[i+1]-_x[i])*(_x[i+1]-_x[i]) + 12.0*(_y[i+1]+_y[i]))/24.0;
    }
    for(size_t i=0; i<_x.size()-1; i++){
        ls[i+1]+= ls[i];
    }
    for(size_t i=0; i<_x.size(); i++){
        ls[i+1]*=-1;
    }
    return ls;
}

double CubicSpline::sintall(){
    std::vector<double> ls = this->sinta();
    return ls[ls.size()-1];
}

double CubicSpline::spl(const double& xx){
    size_t ixx;
    for(size_t i=0; i<_x.size()-1;i++){
        if(xx > _x[i]-(double)1e-12 && xx < _x[i+1]+(double)1e-12){
            ixx = i;
            break;
        }
    }
    if (xx < _x[0] || xx > _x.back()) {
        throw std::out_of_range("xx is outside the spline domain");
    }
    double Ci = _M[ixx]*(_x[ixx+1]-xx)*(_x[ixx+1]-xx)*(_x[ixx+1]-xx)/(double)6.0/_h[ixx] - 
        _M[ixx+1]*(_x[ixx]-xx)*(_x[ixx]-xx)*(_x[ixx]-xx)/(double)6.0/_h[ixx] + 
        (_y[ixx] - _M[ixx]*_h[ixx]*_h[ixx]/(double)6.0)*(_x[ixx+1]-xx)/_h[ixx] 
        + (_y[ixx+1] - _M[ixx+1]*_h[ixx]*_h[ixx]/(double)6.0)*(xx-_x[ixx])/_h[ixx];
    return Ci;
}
