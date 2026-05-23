#include "../include/SchemerBuilder.h"
#include <iostream>

SchemerBuilder& SchemerBuilder::set_alpha(const double value) {
    physics.alpha = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_beta(const double value) {
    physics.beta = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_sigma(const double value) {
    physics.sigma = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_dt(const double value) {
    solving.dt = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_grid_size(const long long int value) {
    solving.grid_points = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_length(const double value) {
    physics.length = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_drift(const double value) {
    if (force_type != ForceType::Drift)
        std::cerr << "[WARNING]: force already set, drift will be used instead\n";
    physics.force_mode = value;
    force_type = ForceType::Drift;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_theta(const double value) {
    solving.theta = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_verbosity(const int value) {
    verbose = value;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_initial_conditions(const Eigen::VectorXd &value) {
    initial_vector = value;
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_initial_conditions(const std::vector<double> &value) {
    initial_vector = Eigen::VectorXd::Map(value.data(), static_cast<Eigen::Index>(value.size()));
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_zero_point(const long long int value) {
    zero_index = value;
    initialization_type = InitialType::Dirac;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(const std::vector<double> &value) {
    if (force_type == ForceType::Drift && std::get<double>(physics.force_mode) != 0.0)
        std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
    force_vector = value;
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(std::function<double(double)> value) {
    if (force_type == ForceType::Drift && std::get<double>(physics.force_mode) != 0.0)
        std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
    force_function = std::move(value);
    force_type = ForceType::Function;
    return *this;
}

Schemer SchemerBuilder::build() {
    if (solving.dt <= 0.0)
        throw std::invalid_argument("dt must be positive");

    if (solving.grid_points <= 0)
        throw std::invalid_argument("I must be positive");

    if (physics.length <= 0.0)
        throw std::invalid_argument("length must be positive");

    if (physics.alpha <= 0.0 || physics.alpha > 2.0)
        throw std::invalid_argument("alpha must be in (0, 2]");

    if (solving.theta < 0.0 || solving.theta > 1.0)
        throw std::invalid_argument("theta must be in [0, 1]");

    if (std::abs(physics.beta) > 1.0)
        throw std::invalid_argument("beta must be in [-1, 1]");

    const double dx = physics.length / solving.grid_points;
    Eigen::Index starting_index;

    Eigen::VectorXd starting_values;
    switch (initialization_type) {
        case InitialType::Middle: {
            starting_index = solving.grid_points / 2;
            starting_values = Eigen::VectorXd::Zero(solving.grid_points+1);
            starting_values(starting_index) = static_cast<double>(solving.grid_points);
            break;
        }
        case InitialType::Dirac: {
            if (zero_index < 0 || zero_index > solving.grid_points)
                throw std::invalid_argument("initial_index must be in [0, I]");
            starting_index = zero_index;
            starting_values = Eigen::VectorXd::Zero(solving.grid_points+1);
            starting_values(zero_index) = static_cast<double>(solving.grid_points);
            break;
        }
        case InitialType::Vector: {
            if (initial_vector.size() != solving.grid_points + 1)
                throw std::invalid_argument("initial values must have size I+1");
            if (zero_index < 0 || zero_index > solving.grid_points)
                throw std::invalid_argument("initial_index must be in [0, I]");
            // if (zero_index < 0)
            //     starting_index = solving.grid_points/ 2;
            starting_index = zero_index;
            starting_values = initial_vector;
            break;
        }
        default:
            throw std::logic_error("Unhandled InitialType enum value");
    }

    switch (force_type) {
        case ForceType::Drift:
            // force_values.assign(grid_size + 1, 0.0);
            break;
        case ForceType::Function: {
            if (force_function == nullptr)
                throw std::invalid_argument("force function not set");
            std::vector<double> force_values;
            force_values.reserve(solving.grid_points + 1);
            for (Eigen::Index i = 0; i <= solving.grid_points; i++) {
                const double position = dx * (i - starting_index);
                force_values.push_back(force_function(position));
            }
            physics.force_mode = std::move(force_values);
            break;
        }
        case ForceType::Vector: {
            if (force_vector.size() != solving.grid_points + 1)
                throw std::invalid_argument("force vector must have size I+1");
            physics.force_mode = std::move(force_vector);
            break;
        }
        default:
            throw std::logic_error("Unhandled ForceType enum value");
    }


    return {physics, solving, starting_values, verbose};
}