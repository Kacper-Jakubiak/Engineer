#include "../include/SchemerBuilder.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "../include/FractionalScheme.h"

void SchemerBuilder::validate_parameters() const {
    if (params.alpha <= 0.0 || params.alpha > 2.0) throw std::invalid_argument("alpha must be in (0, 2]");
    if (std::abs(params.beta) > 1.0) throw std::invalid_argument("beta must be in [-1, 1]");
    if (params.sigma < 0.0) throw std::invalid_argument("sigma must be non-negative");
    if (params.length <= 0.0) throw std::invalid_argument("length must be positive");

    if (params.dt <= 0.0) throw std::invalid_argument("dt must be positive");
    if (params.num_intervals <= 0) throw std::invalid_argument("grid_points must be positive");
    if (params.theta < 0.0 || params.theta > 1.0) throw std::invalid_argument("theta must be in [0, 1]");
    if (params.log_interval_percent < 0.0) throw std::invalid_argument("log interval percent must be non-negative");

    // (beta != 0) is not supported when alpha = 1.0
    if (std::abs(params.alpha - 1) < FractionalScheme::ALPHA_EPSILON && params.beta != 0.0) {
        throw std::invalid_argument("beta != 0.0 not supported for alpha = 1.0");
    }

    // Soft warning for large time steps
    if (params.dt > 1.0) {
        std::cerr << "[WARNING] dt > 1.0\n" << std::endl;
    }
}

Eigen::Index SchemerBuilder::compute_starting_index() const {
    switch (zero_type) {
        case ZeroType::Middle:
            // Place zero reference point at the exact center of the grid
            return params.num_intervals / 2;

        case ZeroType::Index:
            // Use the explicit user-provided grid index
            if (zero_index < 0 || zero_index > params.num_intervals)
                throw std::invalid_argument("initial_index must be in [0, I]");
            return zero_index;

        case ZeroType::Distance: {
            // Convert physical distance offset to nearest discrete grid index
            if (zero_distance < 0.0 || zero_distance > params.length)
                throw std::invalid_argument("initial_distance must be in [0, length]");

            const auto grid_points_double = static_cast<double>(params.num_intervals);
            const double ideal_index = grid_points_double * zero_distance / params.length;
            const auto real_index = static_cast<Eigen::Index>(std::round(ideal_index));
            const double actual_distance = static_cast<double>(real_index) * params.length / grid_points_double;
            const double distance_difference = zero_distance - actual_distance;

            // Warn if distance doesn't align perfectly with grid points
            if (std::abs(distance_difference) > 1e-9) {
                std::cerr << "[WARNING] Target distance " << zero_distance
                        << " snapped to grid index " << real_index
                        << " (discrepancy: " << distance_difference << ")\n";
            }

            if (real_index < 0 || real_index > params.num_intervals)
                throw std::logic_error("Calculated grid index is outside the grid boundaries");

            return real_index;
        }

        default:
            throw std::logic_error("Unhandled ZeroType enum value");
    }
}

Eigen::VectorXd SchemerBuilder::compute_initial_state(const Eigen::Index starting_index) const {
    Eigen::VectorXd starting_values = Eigen::VectorXd::Zero(params.num_intervals + 1);

    switch (initialization_type) {
        case InitialType::Dirac:
            // Scale Dirac delta point mass based on spatial resolution (1 / dx)
            starting_values(starting_index) = static_cast<double>(params.num_intervals) / params.length;
            break;

        case InitialType::Vector:
            // Copy custom user-defined initial vector directly
            if (initial_vector.size() != params.num_intervals + 1)
                throw std::invalid_argument("initial values must have size I+1");
            starting_values = initial_vector;
            break;

        default:
            throw std::logic_error("Unhandled InitialType enum value");
    }

    return starting_values;
}

std::vector<double> SchemerBuilder::compute_force_values(const Grid& grid) const {
    std::vector<double> force_values;

    switch (force_type) {
        case ForceType::Drift:
            // Constant drift uses double in force_mode variant, no vector needed
            break;

        case ForceType::Function: {
            // Evaluate continuous force function at each discrete grid point
            if (!force_function)
                throw std::invalid_argument("force function not set");

            force_values.resize(grid.size);
            for (Eigen::Index i = 0; i < grid.size; ++i) {
                const double position = grid.coordinates[i];
                force_values[i] = force_function(position);
            }
            break;
        }

        case ForceType::Vector:
            // Use pre-computed spatial force vector
            if (force_vector.size() != grid.size)
                throw std::invalid_argument("force vector must have size I+1");
            force_values = force_vector;
            break;

        default:
            throw std::logic_error("Unhandled ForceType enum value");
    }
    return force_values;
}

Grid SchemerBuilder::build_grid(const Eigen::Index starting_index, const double dx) const {
    Grid grid;
    grid.size = params.num_intervals + 1;
    grid.dx = dx;
    grid.starting_index = starting_index;

    grid.coordinates.resize(grid.size);
    for (Eigen::Index i = 0; i < grid.size; ++i) {
        grid.coordinates[i] = dx * static_cast<double>(i - starting_index);
    }
    return grid;
}

Eigen::MatrixXd SchemerBuilder::assemble_step_matrix(const Params &local_params, const Grid &grid) const {
    const FractionalScheme fractional_scheme{local_params, grid};

    Eigen::MatrixXd diffusion_matrix = fractional_scheme.build_diffusion_matrix();

    Eigen::MatrixXd force_matrix = fractional_scheme.build_force_matrix();

    const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(grid.size, grid.size);

    const Eigen::MatrixXd Lhs = id + params.theta * (force_matrix + diffusion_matrix);
    const Eigen::MatrixXd Rhs = id - (1 - params.theta) * (force_matrix + diffusion_matrix);

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    return solver.solve(Rhs);
}

PreparedSystem SchemerBuilder::build_system() const {
    // 1. Validate all inputs
    validate_parameters();

    Params localized_params = this->params;

    // 2. Resolve spatial variables
    const Eigen::Index starting_index = compute_starting_index();
    const double dx = params.length / static_cast<double>(params.num_intervals);

    Grid grid = build_grid(starting_index, dx);

    if (force_type != ForceType::Drift) {
        localized_params.force_mode = compute_force_values(grid);
    }

    Eigen::MatrixXd step_matrix = assemble_step_matrix(localized_params, grid);

    // 6. Return Immutable Prepared System
    return PreparedSystem{
        std::move(step_matrix),                     // The assembled M matrix
        compute_initial_state(starting_index),      // The starting vector u0
        std::move(grid),                            // The spatial mesh
        std::move(localized_params)                 // The finalized settings
    };
}

SchemerRunner SchemerBuilder::build() const {
    return SchemerRunner(build_system());
}