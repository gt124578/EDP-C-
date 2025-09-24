#pragma once
#include <vector>
#include "Mesh1D.hpp"

void solve_heat_equation(
    std::vector<std::vector<double>>& u, 
    const Mesh1D& mesh, 
    int nt, 
    double dt
);