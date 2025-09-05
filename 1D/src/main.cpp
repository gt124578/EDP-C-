#include "parameters.hpp"
#include "solver.hpp"
#include <vector>
#include <cmath>
#include <fstream>

int main() {
    std::vector<std::vector<double>> u(nt, std::vector<double>(nx, 0.0));

    // Condition initiale
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        u[0][i] = std::sin(M_PI * x);
    }

    // Résolution
    solve_heat_equation(u);

    // Sauvegarde dans un CSV pour tracer ensuite
    std::ofstream file("data/solution.csv");
    for (int n = 0; n < nt; ++n) {
        for (int i = 0; i < nx; ++i) {
            file << u[n][i];
            if (i < nx - 1) file << ",";
        }
        file << "\n";
    }
    file.close();

    return 0;
}
