#pragma once
#include "MatrixLUT.h"

class Solver1d {

public:
    Solver1d(Index I, double L, double d, Index J, double t, double alpha);
    void solve();
private:
    MatrixLUT mLUT;
    Index I;
    Index J;
    double L;
    double d;
    double t;
};