#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
Index I = 1000;

double get_lambda(Index k, Index m, double gamma) {
    if (k < 0 || m <= 0)
        throw std::runtime_error("Invalid lambda invocation");


    // if (k == 0 || k == I) {
    //     return pow(m-1, gamma+1) - (m - gamma - 1) * pow(m, gamma);
    // }

    return pow(m+1, gamma+1) - 2*pow(m, gamma+1) + pow(m-1, gamma+1);
}


MatrixXd calculate_V(bool sign) {
    MatrixXd V = MatrixXd::Zero(I+1,I+1);
    for (Index i = 0; i < I+1; i++) {
        V(i, i) = 3.0;
        if (i >= 1)
            V(i, i-1) = -4.0;
        if (i >= 2)
            V(i, i-2) = 1.0;
    }
    // if (!sign)
    //     return -V.transpose();
    return V;
}


MatrixXd calculate_R1(double gamma) {
    MatrixXd R = MatrixXd::Zero(I+1,I+1);
    for (Index i = 0; i < I+1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, i-k, gamma);
            R(i, k) += 3.0 * lambda;
            if (k >= 1)
                R(i, k-1) += -4.0 * lambda;
            if (k >= 2)
                R(i, k-2) += 1.0 * lambda;
        }
    }
    return -1.0/2 * R;
}

MatrixXd calculate_R1t(double gamma) {
    MatrixXd R = MatrixXd::Zero(I+1,I+1);
    for (Index i = 0; i < I+1; i++) {
        for (Index k = i; k <= I; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, k-i, gamma);
            R(i, k) += -3.0 * lambda;
            if (k + 1 <= I)
                R(i, k+1) += 4.0 * lambda;
            if (k + 2 <= I)
                R(i, k+2) += -1.0 * lambda;
        }
    }
    return -1.0/2 * R;
}

int main() {
    IOFormat fmt(
    2,
    DontAlignCols,
    ",",
    "\n",
    "[", "]"
);
    cout << std::fixed;

    double alpha = 0.75;
    constexpr double pi = std::numbers::pi;
    double length = 20.0;
    double dt = 0.01;
    double K = 1.0;
    double mi = 1.0;
    double theta = 0.0;
    int J = 10;

    double gamma = std::ceil(alpha) - alpha;
    double dx = length / static_cast<double>(I);

    const double denominator = 2.0 * std::cos(pi * alpha / 2.0) * std::tgamma(2.0 + gamma) *
                                   std::pow(dx, alpha);
    double omega = -K * dt / denominator;
    double ni = mi * dt / (2.0 * dx);

    VectorXd f = VectorXd::Zero(I + 1);
    f(I / 2) = length / dx;
    const MatrixXd Id = MatrixXd::Identity(I + 1, I + 1);
    const MatrixXd V = calculate_V(true);
    const MatrixXd R = calculate_R1(gamma);
    const MatrixXd Rt = calculate_R1t(gamma);

    // cout << V.format(fmt) << endl << endl;

    const MatrixXd A = R - Rt;
    const MatrixXd lhs = Id + theta * (ni*V + omega * A);
    const MatrixXd rhs = Id - (1.0 - theta) * (ni*V + omega * A);
    // cout << V << endl;
    // cout << A << endl;


    LDLT<MatrixXd> solver(lhs);
    for (int h = 0; h < J; h++) {
        // cout << f.transpose().format(fmt) << endl;
        cout << f.sum() << endl;
        const VectorXd next = solver.solve(rhs * f);
        f = next;
    }
    std::cout << f.transpose().format(fmt) << std::endl;
}