#include "Schemer.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>

Schemer::Schemer(double alpha, double beta, double K, double dt, Eigen::Index I, double length, double mi,
                 double theta, int verbose) : alpha(alpha), beta(beta), K(K), dt(dt), I(I), length(length),
                                              mi(mi), theta(theta), verbose(verbose) {
    dx = length / static_cast<double>(I);
    n = std::ceil(alpha);
    ni = mi * dt / (2.0 * dx);

    if (std::abs(alpha - 1.0) < 1e-7)
        state = 3;
    else if (alpha < 1.0)
        state = 1;
    else
        state = 2;

    if (state == 3) {
        omega = K * dt / (dx * std::numbers::pi);
        L = 0.0;
        R = 0.0;
    } else {
        const double denominator = 2.0 * std::cos(std::numbers::pi * alpha / 2.0);
        L = -(1.0 + beta) / denominator;
        R = -(1.0 - beta) / denominator;
        omega = K * dt / (std::tgamma(n + 1 - alpha) * pow(dx, alpha));
    }

    if (verbose > 0) {
        std::cout << "state " << state << std::endl;
    }
}


double Schemer::lambda(Eigen::Index diagonal) const {
    return pow(diagonal, n - alpha) - pow(diagonal - 1, n - alpha);
}


Eigen::MatrixXd Schemer::get_diffusion_M1() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(I + 1, I + 1);

    double v0 = lambda(1);
    M.diagonal().setConstant(-1.0 * (L + R) * omega * v0);

    for (Eigen::Index i = 1; i <= I; i++) {
        double value = lambda(i) - lambda(i + 1);
        value *= omega;
        M.diagonal(-i).setConstant(value * L);
        M.diagonal(i).setConstant(value * R);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_diffusion_M2() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(I + 1, I + 1);

    double v0 = lambda(2) - 2.0 * lambda(1);
    M.diagonal().setConstant(-1.0 * (L + R) * omega * v0);

    double v1 = lambda(3) - 2 * lambda(2) + lambda(1);
    double v2 = lambda(1);
    M.diagonal(-1).setConstant(-1.0 * (L * omega * v1 + R * omega * v2));
    M.diagonal(1).setConstant(-1.0 * (R * omega * v1 + L * omega * v2));

    for (Eigen::Index i = 2; i <= I; i++) {
        double value = lambda(i + 2) - 2 * lambda(i + 1) + lambda(i);
        value *= omega;
        M.diagonal(-i).setConstant(-1.0 * value * L);
        M.diagonal(i).setConstant(-1.0 * value * R);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_diffusion_M3() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(I + 1, I + 1);

    M.diagonal().setConstant(2.0 * omega);

    for (Eigen::Index i = 1; i <= I; i++) {
        double value = 1.0 / (2.0 * i + 1.0) - 1.0 / (2.0 * (i - 1) + 1.0);
        value *= omega;
        M.diagonal(-i).setConstant(value);
        M.diagonal(i).setConstant(value);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_drift() const {
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


void Schemer::run(Eigen::Index J) {
    using namespace Eigen;
    const MatrixXd id = MatrixXd::Identity(I + 1, I + 1);
    MatrixXd diffusion;

    switch (state) {
        case 1: diffusion = get_diffusion_M1();
            break;
        case 2: diffusion = get_diffusion_M2();
            break;
        case 3: diffusion = get_diffusion_M3();
            break;
        default: std::cout << "Wrong state: " << state << std::endl;
            exit(1);
    }

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


void Schemer::save_history(const std::string &filename) const {
    Eigen::IOFormat fmt(4, Eigen::DontAlignCols, ";", "\n", "");
    const std::string path = std::string(PROJECT_ROOT) + "/" + filename;
    std::ofstream out(path);
    for (const Eigen::VectorXd &v: history) {
        out << v.transpose().format(fmt) << std::endl;
    }
    std::cout << "Saved to " << path << std::endl;
}
