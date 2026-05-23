#include "../include/Schemer.h"

void Schemer::initialize_params() {
    const double alpha = physics.alpha;

    dx = physics.length / static_cast<double>(solving.grid_points);
    n = std::ceil(alpha);
    size = solving.grid_points + 1;

    if (std::holds_alternative<double>(physics.force_mode)) {
        const double mi_temp = std::get<double>(physics.force_mode);
        mi = mi_temp;
        force_type = ForceType::Drift;
    }
    else {
        auto& force_temp = std::get<std::vector<double>>(physics.force_mode);
        force = std::move(force_temp);
        force_type = ForceType::Force;
    }


    if (std::abs(alpha - 1.0) < 1e-7)
        state = SchemeType::Cauchy;
    else if (alpha < 1.0)
        state = SchemeType::Subdiffusive;
    else
        state = SchemeType::Superdiffusive;

    if (state == SchemeType::Cauchy) {
        omega = physics.sigma * solving.dt / (dx * std::numbers::pi);
        L = 0.0;
        R = 0.0;
    } else {
        const double denominator = 2.0 * std::cos(std::numbers::pi * alpha / 2.0);
        L = -(1.0 + physics.beta) / denominator;
        R = -(1.0 - physics.beta) / denominator;
        omega = physics.sigma * solving.dt / (std::tgamma(n + 1 - alpha) * pow(dx, alpha));
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
    }

    const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(size, size);

    const Eigen::MatrixXd Lhs = id + solving.theta * (force_matrix + diffusion_matrix);
    const Eigen::MatrixXd Rhs = id - (1 - solving.theta) * (force_matrix + diffusion_matrix);

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    step_matrix = solver.solve(Rhs);
}
