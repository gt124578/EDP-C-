#include "boundary.hpp"
using namespace std;

//1D
void apply_boundary_conditions(std::vector<double>& u_n) {
    u_n.front()=0.0;
    u_n.back()=0.0;
}

//2D
void apply_boundary_conditions_2D(std::vector<std::vector<double>>& u_n) {
    int nx=u_n.size();
    if (nx==0) return;
    int ny=u_n[0].size();

    //Bords gauche et droit (pour chaque ligne y)
    for (int j=0; j <ny; ++j) {
        u_n[0][j]=0.0;
        u_n[nx-1][j]=0.0;
    }
    //Bords haut et bas (pour chaque colonne x)
    for (int i=0; i<nx;++i) {
        u_n[i][0]=0.0;
        u_n[i][ny - 1]=0.0;
    }
}