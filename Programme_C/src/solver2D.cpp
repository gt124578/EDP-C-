#include "solver2D.hpp"
#include "parameters.hpp"
#include "boundary.hpp"
#include "source_term.hpp"
#include <cmath>

void Solver2D::solve(SolutionHistory& u,const Mesh2D& mesh,int nt,double dt) {
    const int nx=mesh.size_x();
    const int ny=mesh.size_y();
    const double dx=mesh.dx();
    const double dy=mesh.dy();

    //Boucle temporelle
    for (int n=0;n<nt-1;++n) {
        for (int i=1;i<nx-1;++i) {
            for (int j=1;j<ny-1;++j) {
                //Diffusion
                double diffusion_x=(u[n][i+1][j]-2*u[n][i][j]+u[n][i-1][j])/(dx*dx);
                double diffusion_y=(u[n][i][j+1]-2*u[n][i][j]+u[n][i][j-1])/(dy*dy);

                //Convection
                double convection_x=(u[n][i+1][j]-u[n][i-1][j])/(2*dx);
                double convection_y=(u[n][i][j+1]-u[n][i][j-1])/(2*dy);

                //Terme source
                double source=source_term_2D(mesh.x()[i],mesh.y()[j],n*dt);

                //Mise à jour
                u[n+1][i][j]=u[n][i][j]+dt*(
                    D_COEFF*(diffusion_x+diffusion_y)
                    -C_COEFF_X*convection_x
                    -C_COEFF_Y*convection_y
                    +source
                );
            }
        }
        apply_boundary_conditions_2D(u[n+1]);
    }
}