#pragma once

#include "Schemer.h"
#include <functional>
#include <vector>
#include <variant>

class SchemerBuilder {
private:
    constexpr static std::string_view PARAMETER_LOG_FILE_PATH = "parameters.log";

    enum class ZeroType {
        Middle,
        Index,
        Distance
    };

    enum class InitialType {
        Dirac,
        Vector,
    };

    enum class ForceType {
        Drift,
        Function,
        Vector
    };

    Params params = {
        .alpha = 0.0,
        .beta = 0.0,
        .sigma = 1.0,
        .length = 40.0,
        .force_mode = 0.0,

        .dt = 0.001,
        .num_intervals = 1000,
        .theta = 0.5,

        .delimiter = "\t",
        .log_interval_percent = 5.0,
        .verbose = 0,

    };

    InitialType initialization_type = InitialType::Dirac;
    Eigen::VectorXd initial_vector;

    ZeroType zero_type = ZeroType::Middle;
    Eigen::Index zero_index = -1;
    double zero_distance = 0.0;

    ForceType force_type = ForceType::Drift;
    std::function<double(double)> force_function;
    std::vector<double> force_vector;

    void validate_parameters() const;

    [[nodiscard]] Eigen::Index compute_starting_index() const;

    [[nodiscard]] Eigen::VectorXd compute_initial_state(Eigen::Index starting_index) const;

    [[nodiscard]] std::vector<double> compute_force_values(Eigen::Index starting_index, double dx) const;

public:
    SchemerBuilder &set_alpha(double value);

    SchemerBuilder &set_beta(double value);

    SchemerBuilder &set_sigma(double value);

    SchemerBuilder &set_length(double value);

    SchemerBuilder &set_drift(double value);

    SchemerBuilder &set_force(std::vector<double> value);

    SchemerBuilder &set_force(std::function<double(double)> value);

    SchemerBuilder &set_dt(double value);

    SchemerBuilder &set_num_intervals(Eigen::Index value);

    SchemerBuilder &set_theta(double value);

    SchemerBuilder &set_verbosity(int value);

    SchemerBuilder &set_delimiter(std::string value);

    SchemerBuilder &set_log_interval_percent(int value);

    SchemerBuilder &set_initial_conditions(const Eigen::VectorXd &value);

    SchemerBuilder &set_initial_conditions(const std::vector<double> &value);

    SchemerBuilder &set_zero_index(Eigen::Index value);

    SchemerBuilder &set_zero_distance(double value);

    SchemerBuilder &set_zero_middle();

    [[nodiscard]] Schemer build() const;

    [[nodiscard]] Params build_params() const;
};
