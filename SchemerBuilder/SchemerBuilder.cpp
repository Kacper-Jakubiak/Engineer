#include "../include/SchemerBuilder.h"

SchemerBuilder& SchemerBuilder::set_alpha(const double value) {
    alpha = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_beta(const double value) {
    beta = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_K(const double value) {
    K = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_dt(const double value) {
    dt = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_I(const Eigen::Index value) {
    I = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_length(const double value) {
    length = value;
    return *this;
}

SchemerBuilder& SchemerBuilder::set_mi(const double value) {
    mi = value;
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

SchemerBuilder & SchemerBuilder::set_initial(const Eigen::VectorXd &value) {
    initial_vector = value;
    initialization_type = InitialType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_initial(const Eigen::Index value) {
    dirac_index = value;
    initialization_type = InitialType::Dirac;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(const std::vector<double> &value) {
    force_vector = value;
    force_type = ForceType::Vector;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(double(*value)(double)) {
    force_function = value;
    force_type = ForceType::Function;
    return *this;
}

SchemerBuilder & SchemerBuilder::set_force(const double value) {
    drift = value;
    force_type = ForceType::Drift;
    return *this;
}

Schemer SchemerBuilder::build() const {
    if (dt <= 0.0)
        throw std::invalid_argument("dt must be positive");

    if (I <= 0)
        throw std::invalid_argument("I must be positive");

    if (length <= 0.0)
        throw std::invalid_argument("length must be positive");

    if (alpha <= 0.0 || alpha > 2.0)
        throw std::invalid_argument("alpha must be in (0, 2]");

    if (theta < 0.0 || theta > 1.0)
        throw std::invalid_argument("theta must be in [0, 1]");

    if (std::abs(beta) > 1.0)
        throw std::invalid_argument("beta must be in [-1, 1]");

    const double dx = length / I;

    Eigen::VectorXd starting_values;
    switch (initialization_type) {
        case InitialType::Middle:
            starting_values = Eigen::VectorXd::Zero(I+1);
            starting_values(I/2) = static_cast<double>(I);
            break;
        case InitialType::Dirac:
            if (dirac_index < 0 || dirac_index > I)
                throw std::invalid_argument("initial_index must be in [0, I]");
            starting_values = Eigen::VectorXd::Zero(I+1);
            starting_values(dirac_index) = static_cast<double>(I);
            break;
        case InitialType::Vector:
            if (initial_vector.size() != I + 1)
                throw std::invalid_argument("initial values must have size I+1");
            starting_values = initial_vector;
            break;
    }

    std::vector<double> force_values(I+1, 0.0);
    switch (force_type) {
        case ForceType::Drift:
            force_values.assign(I+1, drift);
            break;
        case ForceType::Function:
            if (force_function == nullptr)
                throw std::invalid_argument("force function not set");
            for (Eigen::Index i = 0; i <= I; i++) {
                const double position = dx * i;
                force_values[i] = force_function(position);
            }
            break;
        case ForceType::Vector:
            if (force_vector.size() != I + 1)
                throw std::invalid_argument("force vector must have size I+1");
            force_values = force_vector;
            break;
    }

    return {
        alpha,
        beta,
        K,
        dt,
        I,
        length,
        mi,
        theta,
        verbose,
        std::move(starting_values),
        std::move(force_values)
    };
}
