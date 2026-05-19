#include "../include/SchemerBuilder.h"
#include <iostream>

SchemerBuilder& SchemerBuilder::set_alpha(const double value) {
    alpha = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_beta(const double value) {
    beta = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_sigma(const double value) {
    sigma = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_dt(const double value) {
    dt = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_grid_size(const long long int value) {
    grid_size = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_length(const double value) {
    length = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_drift(const double value) {
    if (force_type != ForceType::Drift)
        std::cerr << "[WARNING]: force already set, drift will be used instead\n";
    mi = value;
    force_type = ForceType::Drift;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_theta(const double value) {
    theta = value;
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
    initial_vector = Eigen::Map<const Eigen::VectorXd>(
        value.data(),
        static_cast<Eigen::Index>(value.size())
        );
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_zero_point(const long long int value) {
    zero_index = value;
    initialization_type = InitialType::Dirac;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(const std::vector<double> &value) {
    if (mi != 0.0)
        std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
    force_vector = value;
    mi = 0.0;
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(double(*value)(double)) {
    if (mi != 0.0)
        std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
    force_function = value;
    mi = 0.0;
    force_type = ForceType::Function;
    return *this;
}

Schemer SchemerBuilder::build() const {
    if (dt <= 0.0)
        throw std::invalid_argument("dt must be positive");

    if (grid_size <= 0)
        throw std::invalid_argument("I must be positive");

    if (length <= 0.0)
        throw std::invalid_argument("length must be positive");

    if (alpha <= 0.0 || alpha > 2.0)
        throw std::invalid_argument("alpha must be in (0, 2]");

    if (theta < 0.0 || theta > 1.0)
        throw std::invalid_argument("theta must be in [0, 1]");

    if (std::abs(beta) > 1.0)
        throw std::invalid_argument("beta must be in [-1, 1]");

    const double dx = length / grid_size;
    Eigen::Index starting_index;

    Eigen::VectorXd starting_values;
    switch (initialization_type) {
        case InitialType::Middle:
            starting_index = grid_size / 2;
            starting_values = Eigen::VectorXd::Zero(grid_size+1);
            starting_values(starting_index) = static_cast<double>(grid_size);
            break;
        case InitialType::Dirac:
            if (zero_index < 0 || zero_index > grid_size)
                throw std::invalid_argument("initial_index must be in [0, I]");
            starting_index = zero_index;
            starting_values = Eigen::VectorXd::Zero(grid_size+1);
            starting_values(zero_index) = static_cast<double>(grid_size);
            break;
        case InitialType::Vector:
            if (initial_vector.size() != grid_size + 1)
                throw std::invalid_argument("initial values must have size I+1");
            starting_index = grid_size / 2;
            starting_values = initial_vector;
            break;
        default:
            throw std::logic_error("Unhandled InitialType enum value");
    }

    std::vector<double> force_values;
    switch (force_type) {
        case ForceType::Drift:
            force_values.assign(grid_size + 1, 0.0);
            break;
        case ForceType::Function:
            if (force_function == nullptr)
                throw std::invalid_argument("force function not set");

            force_values.reserve(grid_size + 1);
            for (Eigen::Index i = 0; i <= grid_size; i++) {
                const double position = dx * (i - starting_index);
                force_values.push_back(force_function(position));
            }
            break;
        case ForceType::Vector:
            if (force_vector.size() != grid_size + 1)
                throw std::invalid_argument("force vector must have size I+1");
            force_values = force_vector;
            break;
        default:
            throw std::logic_error("Unhandled ForceType enum value");
    }

    return {
        alpha,
        beta,
        sigma,
        dt,
        grid_size,
        length,
        mi,
        theta,
        verbose,
        std::move(starting_values),
        std::move(force_values)
    };
}
