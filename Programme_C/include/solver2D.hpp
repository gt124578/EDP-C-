#pragma once
#include <vector>
#include "Mesh2D.hpp"

using namespace std;

class Solver2D {
public:
    using SolutionHistory = vector<vector<vector<double>>>;

    void solve(
        SolutionHistory& u,
        const Mesh2D& mesh,
        int nt,
        double dt
    );
};