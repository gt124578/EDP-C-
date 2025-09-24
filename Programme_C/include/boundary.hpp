#pragma once
#include <vector>

using namespace std;
//1D
void apply_boundary_conditions(vector<double>& u_n);
//2D
void apply_boundary_conditions_2D(vector<vector<double>>& u_n);