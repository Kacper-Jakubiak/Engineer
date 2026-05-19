#include "../include/Schemer.h"

void Schemer::initialize_params() {
    dx = length / static_cast<double>(I);
    n = std::ceil(alpha);
    ni = mi * dt / (2.0 * dx);

    if (std::abs(alpha - 1.0) < 1e-7)
        state = SchemeType::Cauchy;
    else if (alpha < 1.0)
        state = SchemeType::Subdiffusive;
    else
        state = SchemeType::Superdiffusive;

    if (state == SchemeType::Cauchy) {
        omega = K * dt / (dx * std::numbers::pi);
        L = 0.0;
        R = 0.0;
    } else {
        const double denominator = 2.0 * std::cos(std::numbers::pi * alpha / 2.0);
        L = -(1.0 + beta) / denominator;
        R = -(1.0 - beta) / denominator;
        omega = K * dt / (std::tgamma(n + 1 - alpha) * pow(dx, alpha));
    }
}

void Schemer::initialize_matrices() {
    Eigen::MatrixXd diffusion_matrix;

    switch (state) {
        case SchemeType::Subdiffusive: diffusion_matrix = get_diffusion_M1();
            break;
        case SchemeType::Superdiffusive: diffusion_matrix = get_diffusion_M2();
            break;
        case SchemeType::Cauchy: diffusion_matrix = get_diffusion_M3();
            break;
        default: throw std::runtime_error("Invalid Scheme Type");
    }

    Eigen::MatrixXd force_matrix;
    if (mi == 0.0) {force_matrix = get_force();}
    else {force_matrix = get_drift();}

    const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(I + 1, I + 1);

    const Eigen::MatrixXd Lhs = id + theta * (force_matrix + diffusion_matrix);
    const Eigen::MatrixXd Rhs = id - (1 - theta) * (force_matrix + diffusion_matrix);

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    step_matrix = solver.solve(Rhs);
}
