#ifndef MESH2D_H
#define MESH2D_H

#include <vector>

using namespace std;

class Mesh2D {
public:
    Mesh2D(double length_x, double length_y, int points_x, int points_y);

    double dx() const;
    double dy() const;
    int size_x() const;
    int size_y() const;
    const vector<double>& x() const;
    const vector<double>& y() const;

private:
    double _Lx, _Ly;
    int _nx, _ny;
    double _dx, _dy;
    vector<double> _x, _y;
};

#endif