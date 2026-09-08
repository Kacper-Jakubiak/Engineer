/**
 * @file SchemerBuilder.h
 * @brief Header file for the SchemerBuilder class. Helps create and set up Schemer objects.
 */

#pragma once

#include "SchemerRunner.h"
#include "Params.h"
#include "Grid1D.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

/**
 * @class SchemerBuilder
 * @brief Builder class to set up and create Schemer objects.
 *
 * Allows setting physical values, grid options, starting conditions,
 * and forces step-by-step before building the solver.
 */
class SchemerBuilder {
private:
    /// File path used to save parameter logs.
    constexpr static std::string_view PARAMETER_LOG_FILE_PATH = "parameters.log";

    /**
     * @enum ZeroType
     * @brief Options for setting where position zero is on the grid.
     */
    enum class ZeroType {
        Middle,   ///< Places zero in the middle of the grid.
        Index,    ///< Uses a specific grid spot number for zero.
        Distance  ///< Places zero at a specific distance from the start.
    };

    /**
     * @enum InitialType
     * @brief Options for setting the starting conditions.
     */
    enum class InitialType {
        Dirac,   ///< Starts with a single point spike (Dirac delta).
        Vector,  ///< Starts using a list of values from the user.
    };

    /**
     * @enum ForceType
     * @brief Options for setting the force applied to the system.
     */
    enum class ForceType {
        Drift,     ///< A constant pushing force.
        Function,  ///< A force calculated by a math function.
        Vector     ///< A list of force values for each grid point.
    };

    Params params; /// Stores the settings.

    InitialType initialization_type = InitialType::Dirac; /// Selected starting condition type.
    Eigen::VectorXd initial_vector;                        /// Custom list of starting values.

    ZeroType zero_type = ZeroType::Middle; /// Selected zero location option.
    Eigen::Index zero_index = -1;           /// Grid spot number for zero location.
    double zero_distance = 0.0;             /// Distance offset for zero location.

    ForceType force_type = ForceType::Drift;        /// Selected force type.
    std::function<double(double)> force_function;   /// Function that calculates force at a position.
    std::vector<double> force_vector;               /// List of force values.

    /**
     * @brief Checks if settings are valid before building.
     * @throws std::invalid_argument or std::runtime_error if settings are wrong or missing.
     */
    void validate_parameters() const;

    /**
     * @brief Finds the grid index where zero is located based on zero_type.
     * @return Eigen::Index The grid spot for zero.
     */
    [[nodiscard]] Eigen::Index compute_starting_index() const;

    /**
     * @brief Creates the starting vector of probabilities.
     * @param starting_index The grid spot for zero position.
     * @return Eigen::VectorXd The starting values.
     */
    [[nodiscard]] Eigen::VectorXd compute_initial_state(Eigen::Index starting_index) const;

    /**
     * @brief Calculates the force at each grid point.
     * @param grid
     * @return std::vector<double> List of calculated force values.
     */
    [[nodiscard]] std::vector<double> compute_force_values(const Grid1D &grid) const;

    [[nodiscard]] Grid1D build_grid(Eigen::Index starting_index, double dx) const;

    [[nodiscard]] Eigen::MatrixXd assemble_step_matrix(const Params &local_params, const Grid1D &grid) const;

public:
    /**
     * @brief Sets the alpha parameter.
     * @param value Alpha value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_alpha(double value);

    /**
     * @brief Sets the beta parameter.
     * @param value Beta value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_beta(double value);

    /**
     * @brief Sets the sigma parameter.
     * @param value Sigma value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_sigma(double value);

    /**
     * @brief Sets the total length of the spatial area.
     * @param value Length value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_length(double value);

    /**
     * @brief Sets a steady force, equal at each point.
     * @param value Drift value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_drift(double value);

    /**
     * @brief Sets the force using a list of values for each grid point.
     * @param value List of force values.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_force(std::vector<double> value);

    /**
     * @brief Sets the force using a function.
     * @param value A function that takes position and returns the force value.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_force(std::function<double(double)> value);

    /**
     * @brief Sets the time step size (dt).
     * @param value Time step size.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_dt(double value);

    /**
     * @brief Sets the number of grid sections.
     * @param value Number of grid sections.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_num_intervals(Eigen::Index value);

    /**
     * @brief Sets the theta parameter for time stepping.
     * @param value Theta value between 0 and 1.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_theta(double value);

    /**
     * @brief Sets the verbosity level of the returned solver.
     * @param value Detail level (higher means more logs).
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_verbosity(int value);

    /**
     * @brief Sets the character used to separate values in output files.
     * @param value Separator string.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_delimiter(std::string value);

    /**
     * @brief Sets how often to report progress as a percentage.
     * @param value Progress percent step.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_log_interval_percent(double value);

    /**
     * @brief Sets starting values using an Eigen vector.
     * @param value List of starting values.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_initial_conditions(const Eigen::VectorXd &value);

    /**
     * @brief Sets starting values using a standard std::vector.
     * @param value List of starting values.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_initial_conditions(const std::vector<double> &value);

    /**
     * @brief Sets where zero is located using an index.
     * @param value Grid index for zero location.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_zero_index(Eigen::Index value);

    /**
     * @brief Sets where zero is located using a distance offset.
     * @param value Distance for zero location from the start of the grid.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_zero_distance(double value);

    /**
     * @brief Places zero directly in the middle of the domain grid.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_zero_middle();

    /**
     * @brief Sets all parameters at once using an existing Params object.
     * @param value Pre-made Params object.
     * @return SchemerBuilder& Reference to this builder.
     */
    SchemerBuilder &set_params(const Params &value);

    /**
     * @brief Checks options and creates the final Schemer object.
     * @return Schemer Built Schemer solver object.
     */
    [[nodiscard]] SchemerRunner build() const;

    /**
     * @brief Checks options and creates the Params object.
     * @return Params Built Params object.
     */
    [[nodiscard]] PreparedSystem build_system() const;
};