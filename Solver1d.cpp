#include "Solver1d.h"

#include <iostream>

Solver1d::Solver1d(Index I, double L, double d, Index J, double t, double alpha) : mLUT(MatrixLUT(I + 1, alpha)) {
    this->I = I;
    this->L = L;
    this->d = d;
    this->J = J;
    this->t = t;
}

void Solver1d::solve() {
    VectorXd f = VectorXd::Zero(I+1);
    while (true) {
        const VectorXd next = mLUT.A.ldlt().solve(f);
        f = next;
        break;
    }
    std::cout << f << std::endl;
}
