#include <cmath>
#include <functional>
#include <tuple>

struct p{
    double param;
    double u0;
    double beta;
    double Gamma;
    double c1;
    double fi;
};

double func(double beta, struct p params){
    double a = params.param;
    return 1.0/tan(beta) - a + beta;
}

double dfunc(double beta, struct p params){
    return -1.0/sin(beta)/sin(beta) + 1.0;
}

double func1(double gam, struct p params){
    return 2.0*params.u0*cos(gam-params.beta) - params.Gamma*gam/2.0/acos(-1.0) + params.c1 - params.fi;
}

double dfunc1(double gam, struct p params){
    return -2.0*params.u0*sin(gam-params.beta) - params.Gamma/2.0/acos(-1.0);
}

double newton(
    std::function<double(double, struct p)> f,      // Функция f(x, args...)
    std::function<double(double, struct p)> df,     // Производная f'(x, args...)
    double initial_guess,                          // Начальное приближение
    double epsilon,                               // Точность
    size_t max_iterations,                         // Макс. число итераций
    struct p params                                  // Доп. параметры функции
) {
    double x = initial_guess;
    double error = f(x, params);
    size_t iter = 0;

    while (std::abs(error) > epsilon && iter < max_iterations) {
        double derivative = df(x, params);
        if (std::abs(derivative) < 1e-10) {  // Защита от деления на 0
            std::cerr << "Derivative is too small, Newton's method fails.\n";
            break;
        }
        x -= f(x, params) / derivative;
        error = f(x, params);
        if(std::abs(error) > pow(10, 7)){
            std::cerr << "Error is too big.\n";
            break;
        }
        iter++;
    }
    return x;
}

std::vector<double> hilbert(const std::vector<double>& y, const std::vector<double>& x){
    std::vector<double> x2(x.size()-1, 0.0);
    std::vector<double> h(x.size()-1, 0.0);
    std::vector<double> res(x.size(), 0.0);
    for(size_t i=0; i<x.size()-1; i++){
        x2[i] = (x[i] + x[i+1])/2.0;
        h[i] = x[i+1] - x[i];
    }
    auto spl = CubicSpline(x, y);
    spl.assemble();
    std::vector<double> y2(x.size()-1, 0.0);
    for(size_t i=0; i<x.size()-1; i++){
        y2[i] = spl.spl(x2[i]);
    }
    for(size_t i=0; i<x.size(); i++){
        double sum = 0;
        for(size_t j=0; j<x.size()-1; j++){
            sum+= 1.0/2.0/acos(-1.0)*y2[j]/tan((x2[j]-x[i])/2.0)*h[i];
        }
        res[i] = sum;
    }
    return res;
}
