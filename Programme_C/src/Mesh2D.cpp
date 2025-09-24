#include "Mesh2D.hpp"

using namespace std;

Mesh2D::Mesh2D(double length_x, double length_y, int points_x, int points_y)
    : _Lx(length_x), _Ly(length_y), _nx(points_x), _ny(points_y) {
    _dx=_Lx/(_nx-1);
    _dy=_Ly/(_ny-1);

    _x.resize(_nx);
    for (int i=0; i<_nx; ++i) _x[i]=i*_dx;

    _y.resize(_ny);
    for (int j=0; j<_ny; ++j) _y[j]=j*_dy;
}

double Mesh2D::dx() const {return _dx;}
double Mesh2D::dy() const {return _dy;}
int Mesh2D::size_x() const {return _nx;}
int Mesh2D::size_y() const {return _ny;}
const vector<double>& Mesh2D::x() const {return _x;}
const vector<double>& Mesh2D::y() const {return _y;}