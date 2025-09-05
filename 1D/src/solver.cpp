#include "solver.hpp"
#include "parameters.hpp"
#include "boundary.hpp"
#include "source_term.hpp"
#include <cmath>

using namespace std;


void solve_heat_equation(vector<vector<double>>& u){
    for (int n=0;n<nt-1;++n){
        for (int i=0;i<nx-1;++i){
            double convection= -C*(u[n][i+1]-u[n][i-1])/(2.0*dx);
            double diffusion = D*(u[n][i+1]-2*u[n][i]+u[n][i-1])/dx**2;
            double source = source_term(i*dx, n*dt);
            u[n+1][i]=u[n][i]+dt*(convection+diffusion+source);
        }
    apply_boundary_conditions(u[n+1]);
    }
}
