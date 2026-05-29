#include <iostream>

#include "../include/SchemerBuilder.h"

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

SchemerBuilder & SchemerBuilder::set_force(std::vector<double> value) {
    if (force_type == ForceType::Drift && std::holds_alternative<double>(physics.force_mode)) {
        if (std::get<double>(physics.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
        }
    }
    force_vector = std::move(value);
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_force(std::function<double(double)> value) {
    if (force_type == ForceType::Drift && std::holds_alternative<double>(physics.force_mode)) {
        if (std::get<double>(physics.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
        }
    }
    force_function = std::move(value);
    force_type = ForceType::Function;
    return *this;
}