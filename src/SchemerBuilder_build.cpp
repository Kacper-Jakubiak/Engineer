/**
 * @file SchemerBuilder_build.cpp
 * @brief Build and validation logic for the solver.
 */

#include "../include/SchemerBuilder.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "../include/FractionalScheme.h"

void SchemerBuilder::validate_parameters() const {
    if (inner_config.alpha <= 0.0 || inner_config.alpha > 2.0) throw std::invalid_argument("alpha must be in (0, 2]");
    if (std::abs(inner_config.beta) > 1.0) throw std::invalid_argument("beta must be in [-1, 1]");
    if (inner_config.sigma < 0.0) throw std::invalid_argument("sigma must be non-negative");
    if (inner_config.length <= 0.0) throw std::invalid_argument("length must be positive");

    if (inner_config.dt <= 0.0) throw std::invalid_argument("dt must be positive");
    if (inner_config.num_intervals <= 0) throw std::invalid_argument("grid_points must be positive");
    if (inner_config.theta < 0.0 || inner_config.theta > 1.0) throw std::invalid_argument("theta must be in [0, 1]");
    if (inner_config.log_interval_percent < 0.0) throw std::invalid_argument("log interval percent must be non-negative");

    // (beta != 0) is not supported when alpha = 1.0
    if (std::abs(inner_config.alpha - 1) < FractionalScheme::ALPHA_EPSILON && inner_config.beta != 0.0) {
        throw std::invalid_argument("beta != 0.0 not supported for alpha = 1.0");
    }

    // Soft warning for large time steps
    if (inner_config.dt > 1.0) {
        std::cerr << "[WARNING] dt > 1.0\n" << std::endl;
    }
}

Eigen::Index SchemerBuilder::compute_starting_index() const {
    switch (zero_index_type) {
        case ZeroIndexEnum::Middle:
            // Place zero reference point at the exact center of the grid
            return inner_config.num_intervals / 2;

        case ZeroIndexEnum::Index:
            // Use the explicit user-provided grid index
            if (zero_index < 0 || zero_index > inner_config.num_intervals)
                throw std::invalid_argument("initial_index must be in [0, I]");
            return zero_index;

        case ZeroIndexEnum::Distance: {
            // Convert physical distance offset to nearest discrete grid index
            if (zero_distance < 0.0 || zero_distance > inner_config.length)
                throw std::invalid_argument("initial_distance must be in [0, length]");

            const auto grid_points_double = static_cast<double>(inner_config.num_intervals);
            const double ideal_index = grid_points_double * zero_distance / inner_config.length;
            const auto real_index = static_cast<Eigen::Index>(std::round(ideal_index));
            const double actual_distance = static_cast<double>(real_index) * inner_config.length / grid_points_double;
            const double distance_difference = zero_distance - actual_distance;

            // Warn if distance doesn't align perfectly with grid points
            if (std::abs(distance_difference) > 1e-9) {
                std::cerr << "[WARNING] Target distance " << zero_distance
                        << " snapped to grid index " << real_index
                        << " (discrepancy: " << distance_difference << ")\n";
            }

            if (real_index < 0 || real_index > inner_config.num_intervals)
                throw std::logic_error("Calculated grid index is outside the grid boundaries");

            return real_index;
        }

        default:
            throw std::logic_error("Unhandled ZeroIndexEnum enum value");
    }
}

Eigen::VectorXd SchemerBuilder::compute_initial_state(const Eigen::Index starting_index) const {
    const Eigen::Index grid_size = inner_config.get_size();
    Eigen::VectorXd starting_values = Eigen::VectorXd::Zero(grid_size);

    switch (initial_values_type) {
        case InitialValuesEnum::Dirac:
            // Scale Dirac delta point mass based on spatial resolution (1 / dx)
            starting_values(starting_index) = static_cast<double>(inner_config.num_intervals) / inner_config.length;
            break;

        case InitialValuesEnum::Vector:
            // Copy custom user-defined initial vector directly
            if (initial_vector.size() != grid_size)
                throw std::invalid_argument("initial values must have size I+1");
            starting_values = initial_vector;
            break;

        default:
            throw std::logic_error("Unhandled InitialValuesEnum enum value");
    }

    return starting_values;
}

std::variant<double, std::vector<double>> SchemerBuilder::compute_force(const Eigen::Index starting_index) const {
    switch (force_type) {
        case ForceTypeEnum::Drift:
            return drift_value;
        case ForceTypeEnum::Function: {
            if (!force_function)
                throw std::invalid_argument("force function not set");

            const auto coords = inner_config.get_coordinates(starting_index);
            std::vector<double> force_values(coords.size());

            for (size_t i = 0; i < coords.size(); ++i) {
                force_values[i] = force_function(coords[i]);
            }
            return force_values;
        }

        case ForceTypeEnum::Vector:
            // Use pre-computed spatial force vector
            if (force_vector.size() != inner_config.get_size())
                throw std::invalid_argument("force vector must have size I+1");
            return force_vector;

        default:
            throw std::logic_error("Unhandled ForceTypeEnum enum value");
    }
}

Eigen::MatrixXd SchemerBuilder::assemble_step_matrix(const std::variant<double, std::vector<double>>& force_variant) const {
    const FractionalScheme fractional_scheme{inner_config};

    const Eigen::MatrixXd diffusion_matrix = fractional_scheme.build_diffusion_matrix();
    const Eigen::MatrixXd force_matrix = fractional_scheme.build_force_matrix(force_variant);

    const Eigen::Index grid_size = inner_config.get_size();
    const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(grid_size, grid_size);

    const Eigen::MatrixXd Lhs = id + inner_config.theta * (force_matrix + diffusion_matrix);
    const Eigen::MatrixXd Rhs = id - (1 - inner_config.theta) * (force_matrix + diffusion_matrix);

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    return solver.solve(Rhs);
}

SimulationSetup SchemerBuilder::build_setup() const {
    validate_parameters();

    const Eigen::Index starting_index = compute_starting_index();

    auto force_variant = compute_force(starting_index);
    auto initial_state = compute_initial_state(starting_index);

    Eigen::MatrixXd step_matrix = assemble_step_matrix(force_variant);

    return SimulationSetup{
        inner_config,
        std::move(step_matrix),
        std::move(initial_state),
        starting_index,
        std::move(force_variant)
    };
}

SchemerRunner SchemerBuilder::build() const {
    return SchemerRunner(build_setup());
}