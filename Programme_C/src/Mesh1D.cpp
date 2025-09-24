#include "Mesh1D.hpp"

Mesh1D::Mesh1D(double length, int points) : _L(length), _nx(points) {
    _dx=_L/(_nx - 1);
    _x.resize(_nx);
    for (int i=0; i<_nx; ++i) {
        _x[i]=i*_dx;
    }
}

double Mesh1D::dx() const {return _dx;}
int Mesh1D::size() const {return _nx;}
const std::vector<double>& Mesh1D::x() const {return _x;}