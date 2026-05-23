#pragma once

#include "Schemer.h"
#include <functional>

class SchemerBuilder {
private:
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
        .alpha = 0,
        .beta = 0,
        .sigma = 1.0,
        .length = 0.0,
        .force_mode = 0.0
    };

    SolverParams solving = {
        .dt = 0.01,
        .grid_points = 0,
        .theta = 0.5
    };

    InitialType initialization_type = InitialType::Middle;
    Eigen::VectorXd initial_vector;
    long long int zero_index = -1;

    ForceType force_type = ForceType::Drift;
    std::function<double(double)> force_function;
    std::vector<double> force_vector;

    int verbose = 0;

public:
    SchemerBuilder& set_alpha(double value);
    SchemerBuilder& set_beta(double value);
    SchemerBuilder& set_sigma(double value);
    SchemerBuilder& set_dt(double value);
    SchemerBuilder& set_grid_size(long long int value);
    SchemerBuilder& set_length(double value);
    SchemerBuilder& set_drift(double value);
    SchemerBuilder& set_theta(double value);
    SchemerBuilder& set_verbosity(int value);
    SchemerBuilder& set_initial_conditions(const Eigen::VectorXd &value);
    SchemerBuilder& set_initial_conditions(const std::vector<double> &value);
    SchemerBuilder& set_zero_point(long long int value);
    SchemerBuilder& set_force(const std::vector<double> &value);
    SchemerBuilder& set_force(std::function<double(double)> value);

    [[nodiscard]] Schemer build();
};
