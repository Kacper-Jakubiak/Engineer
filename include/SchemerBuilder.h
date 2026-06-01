#pragma once

#include "Schemer.h"
#include <functional>
#include <vector>
#include <variant>

class SchemerBuilder {
private:
    constexpr static std::string_view PARAMETER_LOG_FILE_PATH = "parameters.log";

    enum class InitialType {
        Middle,
        Dirac,
        Vector,
    };
    enum class ForceType {
        Drift,
        Function,
        Vector
    };

    PhysicsParams physics = {
        .alpha = 0.0,
        .beta = 0.0,
        .sigma = 1.0,
        .length = 40.0,
        .force_mode = 0.0
    };

    SolverParams solving = {
        .dt = 0.001,
        .grid_points = 1000,
        .theta = 0.5
    };

    FrontParams front = {
        .delimiter = "\t",
        .log_interval_percent = 5.0,
        .verbose = 0
    };

    InitialType initialization_type = InitialType::Middle;
    Eigen::VectorXd initial_vector;
    Eigen::Index zero_index = -1;

    ForceType force_type = ForceType::Drift;
    std::function<double(double)> force_function;
    std::vector<double> force_vector;

    void save_parameters() const;
    void validate_parameters() const;
    [[nodiscard]] std::pair<Eigen::Index, Eigen::VectorXd> compute_initial_state() const;
    [[nodiscard]] std::vector<double> compute_force_values(Eigen::Index starting_index, double dx) const;

public:
    SchemerBuilder& set_alpha(double value);
    SchemerBuilder& set_beta(double value);
    SchemerBuilder& set_sigma(double value);
    SchemerBuilder& set_dt(double value);
    SchemerBuilder& set_grid_points(Eigen::Index value);
    SchemerBuilder& set_length(double value);
    SchemerBuilder& set_drift(double value);
    SchemerBuilder& set_theta(double value);
    SchemerBuilder& set_verbosity(int value);
    SchemerBuilder& set_delimiter(std::string value);
    SchemerBuilder& set_log_interval_percent(int value);

    SchemerBuilder& set_initial_conditions(const Eigen::VectorXd& value);
    SchemerBuilder& set_initial_conditions(const std::vector<double>& value);
    SchemerBuilder& set_zero_point(Eigen::Index value);
    
    SchemerBuilder& set_force(std::vector<double> value);
    SchemerBuilder& set_force(std::function<double(double)> value);

    [[nodiscard]] Schemer build() const;
};