#include "../include/SchemerBuilder.h"
#include <iostream>

void SchemerBuilder::validate_parameters() const {
    if (solving.dt <= 0.0) throw std::invalid_argument("dt must be positive");

    if (solving.grid_points <= 0) throw std::invalid_argument("grid_points must be positive");

    if (physics.length <= 0.0) throw std::invalid_argument("length must be positive");

    if (physics.alpha <= 0.0 || physics.alpha > 2.0) throw std::invalid_argument("alpha must be in (0, 2]");

    if (solving.theta < 0.0 || solving.theta > 1.0) throw std::invalid_argument("theta must be in [0, 1]");

    if (std::abs(physics.beta) > 1.0) throw std::invalid_argument("beta must be in [-1, 1]");
}


std::pair<Eigen::Index, Eigen::VectorXd> SchemerBuilder::compute_initial_state() const {
    Eigen::Index starting_index = 0;
    Eigen::VectorXd starting_values = Eigen::VectorXd::Zero(solving.grid_points + 1);

    switch (initialization_type) {
        case InitialType::Middle:
            starting_index = solving.grid_points / 2;
            starting_values(starting_index) = static_cast<double>(solving.grid_points);
            break;

        case InitialType::Dirac:
            if (zero_index < 0 || zero_index > solving.grid_points)
                throw std::invalid_argument("initial_index must be in [0, I]");
            starting_index = zero_index;
            starting_values(zero_index) = static_cast<double>(solving.grid_points);
            break;

        case InitialType::Vector:
            if (initial_vector.size() != solving.grid_points + 1)
                throw std::invalid_argument("initial values must have size I+1");
            if (zero_index >= 0 && zero_index <= solving.grid_points)
                starting_index = zero_index;
            else
                starting_index =  (solving.grid_points / 2);
            starting_values = initial_vector;
            break;

        default:
            throw std::logic_error("Unhandled InitialType enum value");
    }
    return {starting_index, starting_values};
}


std::vector<double> SchemerBuilder::compute_force_values(const Eigen::Index starting_index, const double dx) const {
    std::vector<double> force_values;

    switch (force_type) {
        case ForceType::Drift:
            break;

        case ForceType::Function: {
            if (!force_function)
                throw std::invalid_argument("force function not set");

            force_values.reserve(solving.grid_points + 1);
            for (Eigen::Index i = 0; i <= solving.grid_points; ++i) {
                const double position = dx * (i - starting_index);
                force_values.push_back(force_function(position));
            }
            break;
        }
        case ForceType::Vector:
            if (force_vector.size() != solving.grid_points + 1)
                throw std::invalid_argument("force vector must have size I+1");
            force_values = force_vector;
            break;

        default:
            throw std::logic_error("Unhandled ForceType enum value");
    }
    return force_values;
}


Schemer SchemerBuilder::build() const {
    validate_parameters();

    const double dx = physics.length / solving.grid_points;
    auto [starting_index, starting_values] = compute_initial_state();

    PhysicsParams localized_physics = this->physics;

    if (force_type != ForceType::Drift) {
        localized_physics.force_mode = compute_force_values(starting_index, dx);
    }

    return {
        std::move(localized_physics),
        solving,
        std::move(starting_values),
        verbose
    };
}