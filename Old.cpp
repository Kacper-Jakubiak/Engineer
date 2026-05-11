#include "Old.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

using namespace Eigen;
using namespace std;
Index I;

OldSchemer::OldSchemer(double alpha, double K, double dt, Eigen::Index I, double length, double mi,
                 double theta, int verbose) : alpha(alpha), K(K), dt(dt), I(I), length(length),
                                              mi(mi), theta(theta), verbose(verbose) {
    dx = length / static_cast<double>(I);
    ni = mi * dt / (2.0 * dx);
    gamma = std::ceil(alpha) - alpha;

    if (std::abs(alpha - 1.0) < 1e-7)
        state = 3;
    else if (alpha < 1.0)
        state = 1;
    else
        state = 2;

    constexpr double pi = std::numbers::pi;

    double omega;
    if (state != 3) {
        const double denominator = 2.0 * std::cos(pi * alpha / 2.0) * std::tgamma(2.0 + gamma) *
                                   std::pow(dx, alpha);
        omega = -K * dt / denominator;
    } else {
        omega = -K * dt / (pi * dx);
    }

    if (verbose > 0) {
        std::cout << "state " << state << std::endl;
    }
}


double OldSchemer::get_lambda(Index k, Index m) {
    if (k < 0 || m <= 0 || I < k)
        throw std::runtime_error("Invalid lambda invocation");

    // if (k == 0 || k == I) {
    //     return pow(m - 1, gamma + 1) - (m - gamma - 1) * pow(m, gamma);
    // }

    return pow(m + 1, gamma + 1) - 2 * pow(m, gamma + 1) + pow(m - 1, gamma + 1);
}

double OldSchemer::get_lambda(Index m) {
    return 2.0 / (1.0 + static_cast<double>(m));
}

Eigen::MatrixXd OldSchemer::get_drift() const {
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(I + 1, I + 1);
    for (Eigen::Index i = 0; i < I + 1; i++) {
        V(i, i) = 3.0;
        if (i >= 1)
            V(i, i - 1) = -4.0;
        if (i >= 2)
            V(i, i - 2) = 1.0;
    }
    if (mi < 0)
        return -1.0 * ni * V.transpose();
    return ni * V;
}

MatrixXd OldSchemer::get_R2() {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, i - k);
            if (k + 2 <= I)
                R(i, k + 2) += -1.0 * lambda;
            if (k + 1 <= I)
                R(i, k + 1) += 16.0 * lambda;
            R(i, k) += -30.0 * lambda;
            if (k >= 1)
                R(i, k - 1) += 16.0 * lambda;
            if (k >= 2)
                R(i, k - 2) += -1.0 * lambda;
        }
    }
    return -1.0 / 12 * R;
}

MatrixXd OldSchemer::get_R3() {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = get_lambda(i - k);
            R(i, k) += 1.0 * lambda;
            if (k >= 1)
                R(i, k - 1) += -1.0 * lambda;
        }
    }
    return -1.0 * R;
}


MatrixXd OldSchemer::get_R() {
    switch (state) {
        case 1:
              throw std::runtime_error("Somebody stole my R1");
            //return get_R1();
        case 2:
            return get_R2();
        case 3:
            return get_R3();
        default:
            throw std::invalid_argument("Invalid state");
    }
}

void OldSchemer::save_history(const std::string &filename) const {
    Eigen::IOFormat fmt(4, Eigen::DontAlignCols, ";", "\n", "");
    const std::string path = std::string(PROJECT_ROOT) + "/" + filename;
    std::ofstream out(path);
    for (const Eigen::VectorXd &v: history) {
        out << v.transpose().format(fmt) << std::endl;
    }
    std::cout << "Saved to " << path << std::endl;
}

void OldSchemer::run(Eigen::Index J) {
    using namespace Eigen;
    const MatrixXd id = MatrixXd::Identity(I + 1, I + 1);

    const MatrixXd R = get_R();
    const MatrixXd Rt = R.transpose().eval();
    const MatrixXd diffusion = omega * (R + Rt);

    MatrixXd force = get_drift();


    const MatrixXd A = id + theta * (force + diffusion);
    const MatrixXd B = id - (1 - theta) * (force + diffusion);

    VectorXd f = VectorXd::Zero(I + 1);
    f(I / 2) = length / dx;

    history.clear();
    history.reserve(J + 1);
    history.push_back(f);
    std::cout << "SOLVING..." << std::endl;
    const FullPivLU<MatrixXd> solver(A);
    std::cout << "SIMULATING..." << std::endl;
    for (int h = 0; h < J; h++) {
        if (verbose > 2) {
            std::cout << f.transpose() << std::endl;
            std::cout << f.sum() << std::endl;
        }
        const VectorXd next = solver.solve(B * f);
        f = next;
        history.push_back(f);
    }
    if (verbose > 2) std::cout << f.transpose() << std::endl;
}




void run(int argc, char *argv[]) {


    VectorXd f = VectorXd::Zero(I + 1);
    f(I / 2) = length / dx;
    const MatrixXd Id = MatrixXd::Identity(I + 1, I + 1);
    const MatrixXd V = calculate_V(ni);
    const MatrixXd R = get_R(state, gamma);
    const MatrixXd Rt = R.transpose().eval(); //calculate_R1t(gamma);

    const MatrixXd A = R + Rt;
    const MatrixXd lhs = Id + theta * (ni * V + omega * A);
    const MatrixXd rhs = Id - (1.0 - theta) * (ni * V + omega * A);

    vector<VectorXd> history;
    history.push_back(f);
    FullPivLU<MatrixXd> solver(lhs);
    for (int h = 0; h < J; h++) {
        if (verbose > 1) {
            cout << f.transpose().format(fmt) << endl;
            cout << f.sum() << endl;
        }
        const VectorXd next = solver.solve(rhs * f);
        f = next;
        history.push_back(f);
    }
    cout << f.transpose().format(fmt) << endl;
    std::stringstream ss;
    // ss << I << "_" << length << "_" << dt << "_" << K << "_" << mi << "_" << theta << "_" << J << "_" << alpha << ".txt";
    save_history(history, "old");
}