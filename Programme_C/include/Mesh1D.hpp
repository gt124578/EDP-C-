#ifndef MESH1D_H
#define MESH1D_H

#include <vector>

using namespace std;

class Mesh1D {
public:
    Mesh1D(double length, int points);
    double dx() const;
    int size() const;
    const vector<double>& x() const;

private:
    double _L;
    int _nx;
    double _dx;
    vector<double> _x;
};

#endif