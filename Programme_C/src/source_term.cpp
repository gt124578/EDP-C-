#include "source_term.hpp"
#include "parameters.hpp"
#include <cmath>

using namespace std;

//1D
double source_term(double x, double t) {
    return sin(M_PI*x)
         + D_COEFF*(M_PI*M_PI)*sin(M_PI*x)*(1+t) //diffusion
         + C_COEFF*M_PI*cos(M_PI*x)*(1+t); //convection
}

//2D
double source_term_2D(double x,double y,double t) {
    double pi=M_PI;
    double diffusion=2*D_COEFF*pi*pi*sin(pi*x)*sin(pi*y)*(1+t);
    double convection=C_COEFF_X*pi*cos(pi*x)*sin(pi*y)*(1+t)
                    +C_COEFF_Y*pi*sin(pi*x)*cos(pi*y)*(1+t);
    double derivee_t=sin(pi*x)*sin(pi*y);
    return derivee_t+diffusion+convection;
}
