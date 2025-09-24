#include "solver.hpp"
#include "parameters.hpp"
#include "boundary.hpp"
#include "source_term.hpp"
#include <cmath>

using namespace std;


void solve_heat_equation(vector<vector<double>>& u, const Mesh1D& mesh, int nt, double dt) {
    const int nx=mesh.size();
    const double dx=mesh.dx();
    for (int n=0;n<nt-1;++n){
        for (int i=1;i<nx-1;++i){
            double convection=-C_COEFF*(u[n][i+1]-u[n][i-1])/(2.0*dx);
            double diffusion= D_COEFF*(u[n][i+1]-2*u[n][i]+u[n][i-1])/(dx*dx);
            double source=source_term(mesh.x()[i], n*dt);
            u[n+1][i]=u[n][i]+dt*(convection+diffusion+source);
        }
    apply_boundary_conditions(u[n+1]);
    }
}
