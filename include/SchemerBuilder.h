#pragma once

#include "Schemer.h"

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
    double alpha = 0.0;
    double length = 0.0;
    long long int grid_size = 0;

    double beta = 0.0;
    double sigma = 1.0;
    double mi = 0.0;

    double dt = 0.01;
    double theta = 0.5;

    InitialType initialization_type = InitialType::Middle;
    Eigen::VectorXd initial_vector;
    long long int zero_index = 0;

    ForceType force_type = ForceType::Drift;
    double drift = 0.0;
    std::vector<double> force_vector;
    double (*force_function)(double) = nullptr;

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
    SchemerBuilder& set_force(double (*value)(double));

    [[nodiscard]] Schemer build() const;
};
