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
    double alpha = 1.0;
    double beta = 0.0;
    double K = 1.0;
    double dt = 0.01;
    double length = 40.0;
    double mi = 0.0;
    double theta = 0.5;
    Eigen::Index I = 1000;

    InitialType initialization_type = InitialType::Middle;
    Eigen::VectorXd initial_vector;
    Eigen::Index dirac_index = 0;

    ForceType force_type = ForceType::Drift;
    double drift = 0.0;
    std::vector<double> force_vector;
    double (*force_function)(double) = nullptr;

    int verbose = 0;

public:
    SchemerBuilder& set_alpha(double value);
    SchemerBuilder& set_beta(double value);
    SchemerBuilder& set_K(double value);
    SchemerBuilder& set_dt(double value);
    SchemerBuilder& set_I(Eigen::Index value);
    SchemerBuilder& set_length(double value);
    SchemerBuilder& set_mi(double value);
    SchemerBuilder& set_theta(double value);
    SchemerBuilder& set_verbosity(int value);
    SchemerBuilder& set_initial(const Eigen::VectorXd &value);
    SchemerBuilder& set_initial(Eigen::Index value);
    SchemerBuilder& set_force(const std::vector<double> &value);
    SchemerBuilder& set_force(double (*value)(double));
    SchemerBuilder& set_force(double value);

    [[nodiscard]] Schemer build() const;
};
