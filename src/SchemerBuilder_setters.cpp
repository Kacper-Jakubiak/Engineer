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

// --- Force Configuration ---

SchemerBuilder &SchemerBuilder::set_drift(const double value) {
    // Warn the user if a custom force function or vector is being overridden
    if (force_type != ForceType::Drift) {
        std::cerr << "[WARNING]: force was previously set. Drift will be used instead\n";
    }
    params.force_mode = value;
    force_type = ForceType::Drift;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::vector<double> value) {
    // Warn if replacing an existing drift setting
    if (force_type == ForceType::Drift && std::holds_alternative<double>(params.force_mode)) {
        if (std::get<double>(params.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift was previously set. Force will be used instead\n";
        }
    }
    force_vector = std::move(value);
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_force(std::function<double(double)> value) {
    // Warn if replacing an existing drift setting
    if (force_type == ForceType::Drift && std::holds_alternative<double>(params.force_mode)) {
        if (std::get<double>(params.force_mode) != 0.0) {
            std::cerr << "[WARNING]: drift was previously set. Force will be used instead\n";
        }
    }
    force_function = std::move(value);
    force_type = ForceType::Function;
    return *this;
}

// --- Numerical Grid & Stepping ---

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

// --- Logging & Output ---

SchemerBuilder &SchemerBuilder::set_verbosity(const int value) {
    params.verbose = value;
    return *this;
}

SchemerBuilder &SchemerBuilder::set_delimiter(std::string value) {
    params.delimiter = std::move(value);
    return *this;
}

SchemerBuilder &SchemerBuilder::set_log_interval_percent(const double value) {
    params.log_interval_percent = value;
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

// --- Bulk Configuration ---

SchemerBuilder &SchemerBuilder::set_params(const Params &value) {
    // Delegate to existing setter methods to ensure correct state updates
    this->set_alpha(value.alpha);
    this->set_beta(value.beta);
    this->set_sigma(value.sigma);
    this->set_length(value.length);

    // Unpack force variant (drift value vs. force vector)
    if (std::holds_alternative<double>(value.force_mode)) {
        this->set_drift(std::get<double>(value.force_mode));
    } else {
        this->set_force(std::get<std::vector<double>>(value.force_mode));
    }

    this->set_dt(value.dt);
    this->set_num_intervals(value.num_intervals);
    this->set_theta(value.theta);

    this->set_delimiter(value.delimiter);
    this->set_log_interval_percent(value.log_interval_percent);
    this->set_verbosity(value.verbose);

    return *this;
}