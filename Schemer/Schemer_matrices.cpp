#include "../include/Schemer.h"

std::vector<double> Schemer::get_lambdas() const {
    std::vector<double> result(size + 2);
    for (int i = 1; i < size + 2; i++) {
        result[i] = pow(i, n - params.alpha) - pow(i - 1, n - params.alpha);
    }

    return result;
}

Eigen::MatrixXd Schemer::get_diffusion_M1() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(size, size);
    const auto lambda = get_lambdas();

    const double v0 = lambda[1];
    M.diagonal().setConstant(-1.0 * (L + R) * omega * v0);

    for (Eigen::Index i = 1; i < size; i++) {
        double value = lambda[i] - lambda[i + 1];
        value *= omega;
        M.diagonal(-i).setConstant(value * L);
        M.diagonal(i).setConstant(value * R);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_diffusion_M2() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(size, size);
    const auto lambda = get_lambdas();

    double v0 = lambda[2] - 2.0 * lambda[1];
    M.diagonal().setConstant(-1.0 * (L + R) * omega * v0);

    double v1 = lambda[3] - 2 * lambda[2] + lambda[1];
    double v2 = lambda[1];
    M.diagonal(-1).setConstant(-1.0 * (L * omega * v1 + R * omega * v2));
    M.diagonal(1).setConstant(-1.0 * (R * omega * v1 + L * omega * v2));

    for (Eigen::Index i = 2; i < size; i++) {
        double value = lambda[i + 2] - 2 * lambda[i + 1] + lambda[i];
        value *= omega;
        M.diagonal(-i).setConstant(-1.0 * value * L);
        M.diagonal(i).setConstant(-1.0 * value * R);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_diffusion_M3() const {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(size, size);

    M.diagonal().setConstant(2.0 * omega);

    for (Eigen::Index i = 1; i < size; i++) {
        const auto double_i = static_cast<double>(i);
        double value = 1.0 / (2.0 * double_i + 1.0) - 1.0 / (2.0 * (double_i - 1) + 1.0);
        value *= omega;
        M.diagonal(-i).setConstant(value);
        M.diagonal(i).setConstant(value);
    }

    return M;
}


Eigen::MatrixXd Schemer::get_drift() const {
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(size, size);
    for (Eigen::Index i = 0; i < size; i++) {
        V(i, i) = 3.0;
        if (i >= 1)
            V(i, i - 1) = -4.0;
        if (i >= 2)
            V(i, i - 2) = 1.0;
    }

    const double ni = mi * params.dt / (2.0 * dx);

    if (mi < 0)
        return -1.0 * ni * V.transpose();
    return ni * V;
}

Eigen::MatrixXd Schemer::get_force() const {
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(size, size);
    for (Eigen::Index i = 0; i < size; i++) {
        if (i > 0) {
            V(i, i) -= force[i - 1];
            V(i, i - 1) -= force[i];
        }
        if (i < size - 1) {
            V(i, i) += force[i + 1];
            V(i, i + 1) += force[i];
        }
    }
    V *= params.dt / (2 * dx);
    return V;
}