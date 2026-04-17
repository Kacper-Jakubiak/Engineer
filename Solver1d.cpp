#include "Solver1d.h"

#include <iostream>
#include <numbers>
#include <cmath>

Solver1d::Solver1d(double alpha, double theta, Index I, double length, double dt, double K, double mi) :
      theta(theta),
      I(I),
      dx(length / static_cast<double>(I)),
      dt(dt),
      mLUT(I + 1, alpha)
{
    this->ni = mi * dt / (2.0 * dx);

    const double n = std::ceil(alpha);
    constexpr double pi = std::numbers::pi;

    if (std::abs(alpha - 1.0) < 1e-9) {
        this->omega = -K * dt / (pi * dx);
    } else {
        const double denominator = 2.0 * std::cos(pi * alpha / 2.0) * std::tgamma(2.0 + n - alpha) *
                                   std::pow(dx, alpha);
        this->omega = -K * dt / denominator;
    }
}

void Solver1d::solve() const {
    VectorXd f = VectorXd::Zero(I + 1);
    const MatrixXd Id = MatrixXd::Identity(I + 1, I + 1);
    const MatrixXd V = MatrixLUT::calculate_V(I+1, ni);
    while (true) {
        const MatrixXd lhs = Id + theta * (ni*V + omega * mLUT.A);
        const MatrixXd rhs = Id + (1.0 - theta) * (ni*V + omega * mLUT.A);
        const VectorXd next = lhs.ldlt().solve(rhs * f);
        f = next;
        break;
    }
    std::cout << f << std::endl;
}
