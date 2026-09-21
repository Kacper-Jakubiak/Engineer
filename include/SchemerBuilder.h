/**
 * @file SchemerBuilder.h
 * @brief Builder for creating simulation solvers.
 */

#pragma once

#include "SchemerRunner.h"
#include "Config.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

/**
 * @class SchemerBuilder
 * @brief Sets up a simulation step by step.
 *
 * Use this class to configure all settings and create a solver.
 * Each method returns this builder so you can chain calls.
 */
class SchemerBuilder {
private:
    /**
     * @enum ZeroIndexEnum
     * @brief How to set the starting position on the grid.
     */
    enum class ZeroIndexEnum {
        Middle,   ///< Center of grid
        Index,    ///< Specific grid position
        Distance  ///< Position by distance
    };

    /**
     * @enum InitialValuesEnum
     * @brief Type of starting values.
     */
    enum class InitialValuesEnum {
        Dirac,   ///< Single point
        Vector,  ///< Custom values
    };

    /**
     * @enum ForceTypeEnum
     * @brief Type of force applied.
     */
    enum class ForceTypeEnum {
        Drift,     ///< Constant force
        Function,  ///< Force from function
        Vector     ///< Force values
    };

    Config inner_config;

    InitialValuesEnum initial_values_type = InitialValuesEnum::Dirac;
    Eigen::VectorXd initial_vector;

    ZeroIndexEnum zero_index_type = ZeroIndexEnum::Middle;
    Eigen::Index zero_index = -1;
    double zero_distance = 0.0;

    ForceTypeEnum force_type = ForceTypeEnum::Drift;
    double drift_value = 0.0;
    std::function<double(double)> force_function;
    std::vector<double> force_vector;

    /**
     * @brief Check that all settings are valid.
     */
    void validate_parameters() const;

    /**
     * @brief Find the starting grid position.
     * @return Eigen::Index Starting position.
     */
    [[nodiscard]] Eigen::Index compute_starting_index() const;

    /**
     * @brief Create initial values.
     * @param starting_index Starting position.
     * @return Eigen::VectorXd Initial values.
     */
    [[nodiscard]] Eigen::VectorXd compute_initial_state(Eigen::Index starting_index) const;

    /**
     * @brief Compute force values.
     * @param starting_index Starting position.
     * @return Force specification.
     */
    [[nodiscard]] std::variant<double, std::vector<double>> compute_force_variant(Eigen::Index starting_index) const;

    /**
     * @brief Build time-stepping matrix.
     * @param force_variant Force specification.
     * @return Eigen::MatrixXd The matrix.
     */
    [[nodiscard]] Eigen::MatrixXd compute_step_matrix(const std::variant<double, std::vector<double>>& force_variant) const;

public:
    /**
     * @brief Set alpha parameter.
     * @param value Parameter value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_alpha(double value);

    /**
     * @brief Set beta parameter.
     * @param value Parameter value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_beta(double value);

    /**
     * @brief Set sigma parameter.
     * @param value Parameter value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_sigma(double value);

    /**
     * @brief Set domain size.
     * @param value Length value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_length(double value);

    /**
     * @brief Set constant force.
     * @param value Force value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_drift(double value);

    /**
     * @brief Set force from values.
     * @param value Force for each grid point.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_force(std::vector<double> value);

    /**
     * @brief Set force from function.
     * @param value Function for force.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_force(std::function<double(double)> value);

    /**
     * @brief Set time step size.
     * @param value Step size.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_dt(double value);

    /**
     * @brief Set number of grid points.
     * @param value Number of intervals.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_num_intervals(Eigen::Index value);

    /**
     * @brief Set theta parameter.
     * @param value Theta value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_theta(double value);

    /**
     * @brief Set output detail level of the resulting solver.
     * @param value Verbosity level.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_verbosity(int value);

    /**
     * @brief Set output column separator.
     * @param value Separator string.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_delimiter(std::string value);

    /**
     * @brief Set progress log frequency.
     * @param value Percentage (%).
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_log_interval_percent(double value);

    /**
     * @brief Set initial values from Eigen vector.
     * @param value Initial values.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_initial_conditions(const Eigen::VectorXd &value);

    /**
     * @brief Set initial values from std vector.
     * @param value Initial values.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_initial_conditions(const std::vector<double> &value);

    /**
     * @brief Set starting position by grid index.
     * @param value Grid index.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_zero_index(Eigen::Index value);

    /**
     * @brief Set starting position by distance from left side of the grid.
     * @param value Distance value.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_zero_distance(double value);

    /**
     * @brief Set starting position to center.
     * @return SchemerBuilder& This builder.
     */
    SchemerBuilder &set_zero_middle();

    /**
     * @brief Create solver with current settings.
     * @return SchemerRunner Ready solver.
     */
    [[nodiscard]] SchemerRunner build() const;

    /**
     * @brief Get prepared setup.
     * @return RunnerSetup The setup.
     */
    [[nodiscard]] SimulationSetup build_setup() const;
};