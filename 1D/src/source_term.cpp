#include "source_term.hpp"
#include "parameters.hpp"
#include <cmath>

double source_term(double x, double t) {
    return std::sin(M_PI * x)
         + D * (M_PI * M_PI) * std::sin(M_PI * x) * (1 + t) //diffusion
         + C * M_PI * std::cos(M_PI * x) * (1 + t); //convection
}
