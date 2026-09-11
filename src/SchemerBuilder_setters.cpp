/**
 * @file SchemerBuilder_setters.cpp
 * @brief Setter methods for builder configuration.
 */

#include <iostream>
#include <utility>

#include "../include/SchemerBuilder.h"

SchemerBuilder &SchemerBuilder::set_alpha(const double value) {
    config.alpha = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_beta(const double value) {
    config.beta = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_sigma(const double value) {
    config.sigma = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_length(const double value) {
    config.length = value;
    return *this;
}

// --- Force Configuration ---

SchemerBuilder &SchemerBuilder::set_drift(const double value) {
    // Warn the user if a custom force function or vector is being overridden
    if (force_type != ForceType::Drift) {
        std::cerr << "[WARNING]: force was previously set. Drift will be used instead\n";
    }
    drift_value = value;
    force_type = ForceType::Drift;
    force_function = nullptr;
    force_vector.clear();
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::vector<double> value) {
    // Warn if replacing an existing drift setting
    if (force_type == ForceType::Drift && drift_value != 0.0) {
        std::cerr << "[WARNING]: drift was previously set. Force vector will be used instead\n";
    }
    force_vector = std::move(value);
    force_type = ForceType::Vector;
    force_function = nullptr;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::function<double(double)> value) {
    // Warn if replacing an existing drift setting
    if (force_type == ForceType::Drift && drift_value != 0.0) {
        std::cerr << "[WARNING]: drift was previously set. Force function will be used instead\n";
    }
    force_function = std::move(value);
    force_type = ForceType::Function;
    force_vector.clear();
    return *this;
}

// --- Numerical Grid & Stepping ---

SchemerBuilder &SchemerBuilder::set_dt(const double value) {
    config.dt = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_num_intervals(const long long int value) {
    config.num_intervals = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_theta(const double value) {
    config.theta = value;
    return *this;
}

// --- Logging & Output ---

SchemerBuilder &SchemerBuilder::set_verbosity(const int value) {
    config.verbose = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_delimiter(std::string value) {
    config.delimiter = std::move(value);
    return *this;
}

SchemerBuilder &SchemerBuilder::set_log_interval_percent(const double value) {
    config.log_interval_percent = value;
    return *this;
}

// --- Initial Conditions ---

SchemerBuilder &SchemerBuilder::set_initial_conditions(const Eigen::VectorXd &value) {
    initial_vector = value;
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_initial_conditions(const std::vector<double> &value) {
    // Map standard vector data directly into the Eigen vector format
    initial_vector = Eigen::VectorXd::Map(value.data(), static_cast<Eigen::Index>(value.size()));
    initialization_type = InitialType::Vector;
    return *this;
}

// --- Zero Placement Options ---

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