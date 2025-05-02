#include <cmath>
#include <functional>
#include <tuple>

double func(double beta, double param){
    double a = param;
    return 1.0/tan(beta) - a + beta;
}

double dfunc(double beta, double param){
    return -1.0/sin(beta)/sin(beta) + 1.0;
}

double func1(double gam, double beta, double u0, double Gamma, double c1, double fi){
    return 2.8*u0*cos(gam-beta) - Gamma*gam/2/acos(-1.0) + c1 - fi;
}

double dfunc1(double gam, double beta, double u0, double Gamma, double c1, double fi){
    return -2.0*u0*sin(gam-beta) - Gamma/2.0/acos(-1.0);
}

//double newton(double zero, double epsilon, double maxiter, 
        //double param1, double param2=0, double param3=0, double param4=0, double param5=0){
    //double a = param1;
    //double err = func(zero, a);
    //size_t iter = 0;
    //double beta = zero;
    //while(err>epsilon){
        //beta-= func(beta, a)/dfunc(beta, a);
        //if(iter >= maxiter) break;
    //}
    //return beta;
//}

template <typename... Args>
double newton(
    std::function<double(double, Args...)> f,      // Функция f(x, args...)
    std::function<double(double, Args...)> df,     // Производная f'(x, args...)
    double initial_guess,                          // Начальное приближение
    double epsilon,                               // Точность
    size_t max_iterations,                         // Макс. число итераций
    Args... args                                  // Доп. параметры функции
) {
    double x = initial_guess;
    double error = f(x, args...);
    size_t iter = 0;

    while (std::abs(error) > epsilon && iter < max_iterations) {
        double derivative = df(x, args...);
        if (std::abs(derivative) < 1e-10) {  // Защита от деления на 0
            std::cerr << "Derivative is too small, Newton's method fails.\n";
            break;
        }
        x -= f(x, args...) / derivative;
        error = f(x, args...);
        iter++;
    }

    return x;
}
