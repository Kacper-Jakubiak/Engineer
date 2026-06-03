#include "../include/Schemer.h"

void Schemer::initialize_values() {
    const double alpha = params.alpha;

    dx = params.length / static_cast<double>(params.num_intervals);
    n = std::ceil(alpha);
    size = params.num_intervals + 1;

    if (std::holds_alternative<double>(params.force_mode)) {
        force_type = ForceType::Drift;
    }
    else {
        force_type = ForceType::Force;
    }


    if (std::abs(alpha - 1.0) < ALPHA_EPSILON)
        state = SchemeType::Cauchy;
    else if (alpha < 1.0)
        state = SchemeType::Subdiffusive;
    else
        state = SchemeType::Superdiffusive;

    if (state == SchemeType::Cauchy) {
        omega = params.sigma * params.dt / (dx * std::numbers::pi);
        L = 0.0;
        R = 0.0;
    } else {
        const double denominator = 2.0 * std::cos(std::numbers::pi * alpha / 2.0);
        L = -(1.0 + params.beta) / denominator;
        R = -(1.0 - params.beta) / denominator;
        omega = params.sigma * params.dt / (std::tgamma(n + 1 - alpha) * pow(dx, alpha));
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
    switch (force_type) {
        case ForceType::Drift: force_matrix = get_drift();
            break;
        case ForceType::Force: force_matrix = get_force();
            break;
        default:
            throw std::runtime_error("Invalid Force Type");
    }

    const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(size, size);

    const Eigen::MatrixXd Lhs = id + params.theta * (force_matrix + diffusion_matrix);
    const Eigen::MatrixXd Rhs = id - (1 - params.theta) * (force_matrix + diffusion_matrix);

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    step_matrix = solver.solve(Rhs);
}
