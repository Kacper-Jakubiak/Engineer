#include <iostream>
#include <utility>

#include "../include/SchemerBuilder.h"

SchemerBuilder &SchemerBuilder::set_alpha(const double value) {
    params.alpha = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_beta(const double value) {
    params.beta = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_sigma(const double value) {
    params.sigma = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_length(const double value) {
    params.length = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_drift(const double value) {
    if (force_type != ForceType::Drift)
        std::cerr << "[WARNING]: force already set, drift will be used instead\n";
    params.force_mode = value;
    force_type = ForceType::Drift;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::vector<double> value) {
    if (force_type == ForceType::Drift && std::holds_alternative<double>(params.force_mode)) {
        if (std::get<double>(params.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
        }
    }
    force_vector = std::move(value);
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::function<double(double)> value) {
    if (force_type == ForceType::Drift && std::holds_alternative<double>(params.force_mode)) {
        if (std::get<double>(params.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift already set, resetting drift to 0.0. Force will be used instead\n";
        }
    }
    force_function = std::move(value);
    force_type = ForceType::Function;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_dt(const double value) {
    params.dt = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_num_intervals(const long long int value) {
    params.num_intervals = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_theta(const double value) {
    params.theta = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_verbosity(const int value) {
    params.verbose = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_delimiter(std::string value) {
    params.delimiter = std::move(value);
    return *this;
}

SchemerBuilder &SchemerBuilder::set_log_interval_percent(const int value) {
    params.log_interval_percent = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_initial_conditions(const Eigen::VectorXd &value) {
    initial_vector = value;
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_initial_conditions(const std::vector<double> &value) {
    initial_vector = Eigen::VectorXd::Map(value.data(), static_cast<Eigen::Index>(value.size()));
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_zero_index(const long long int value) {
    zero_index = value;
    zero_type = ZeroType::Index;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_zero_distance(const double value) {
    zero_distance = value;
    zero_type = ZeroType::Distance;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_zero_middle() {
    zero_type = ZeroType::Middle;
    return *this;
}