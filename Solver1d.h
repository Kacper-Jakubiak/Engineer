#pragma once
#include "MatrixLUT.h"

class Solver1d {

public:
    Solver1d(double alpha, double theta, Index I, double length, double dt, double K, double mi);
    void solve() const;
private:
    double theta;
    Index I;
    double dx;
    double dt;
    double omega;
    double ni;
    MatrixLUT mLUT;
};