/**
 * @file SchemerBuilder.h
 * @brief Builder class for setting up and creating simulation runners.
 */

#pragma once

#include "SchemerRunner.h"
#include "Config.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <variant>

/**
 * @class SchemerBuilder
 * @brief Configures simulation settings step-by-step and creates a runner.
 *
 * Provides method chaining to set physical parameters, initial grid values,
 * origin position, and external forces before building the final runner object.
 */
class SchemerBuilder {
private:
    /**
     * @enum ZeroIndexEnum
     * @brief Method used to position the spatial origin ($x = 0.0$) on the grid.
     */
    enum class ZeroIndexEnum {
        Middle,   ///< Place origin at the exact center of the grid.
        Index,    ///< Place origin at a specific grid index.
        Distance  ///< Place origin at a physical distance measured from the left edge.
    };

    /**
     * @enum InitialValuesEnum
     * @brief Type of starting values applied across the grid at $t = 0$.
     */
    enum class InitialValuesEnum {
        Dirac,   ///< Single peak (impulse) at the origin ($x = 0.0$).
        Vector   ///< Custom starting values provided for all grid points.
    };

    /**
     * @enum ForceTypeEnum
     * @brief Type of external force applied to the system.
     */
    enum class ForceTypeEnum {
        Drift,     ///< Constant force applied uniformly across the domain.
        Function,  ///< Spatially varying force defined by a function $f(x)$.
        Vector     ///< Custom force value provided for each grid point.
    };

    /// Internal configuration object storing parameters and grid settings.
    Config inner_config;

    /// Selected type for initial starting values.
    InitialValuesEnum initial_values_type = InitialValuesEnum::Dirac;

    /// Custom initial values vector across all grid points.
    Eigen::VectorXd initial_vector;

    /// Selected method for placing origin $x = 0.0$.
    ZeroIndexEnum zero_index_type = ZeroIndexEnum::Middle;

    /// Explicit grid index chosen for $x = 0.0$.
    Eigen::Index zero_index = -1;

    /// Distance from the left edge chosen for $x = 0.0$.
    double zero_distance = 0.0;

    /// Selected force representation type.
    ForceTypeEnum force_type = ForceTypeEnum::Drift;

    /// Constant force value used when force type is set to Drift.
    double drift_value = 0.0;

    /// Mathematical function $f(x)$ used when force varies by position.
    std::function<double(double)> force_function;

    /// List of force values at each grid point.
    std::vector<double> force_vector;

    /**
     * @brief Checks if all settings and parameters are valid before building.
     */
    void validate_parameters() const;

    /**
     * @brief Calculates grid point index corresponding to origin $x = 0.0$.
     * @return Eigen::Index Grid index for the origin.
     */
    [[nodiscard]] Eigen::Index compute_starting_index() const;

    /**
     * @brief Creates starting value vector for all grid points at $t = 0$.
     * @param starting_index Grid point index for origin $x = 0.0$.
     * @return Eigen::VectorXd Initial grid values vector.
     */
    [[nodiscard]] Eigen::VectorXd compute_initial_state(Eigen::Index starting_index) const;

    /**
     * @brief Converts configured external force into scalar or vector format.
     * @param starting_index Grid point index for origin $x = 0.0$.
     * @return std::variant<double, std::vector<double>> Force representation.
     */
    [[nodiscard]] std::variant<double, std::vector<double>> compute_force_variant(Eigen::Index starting_index) const;

    /**
     * @brief Builds time-stepping matrix based on current settings and forces.
     * @param force_variant Computed force representation.
     * @return Eigen::MatrixXd Time step update matrix.
     */
    [[nodiscard]] Eigen::MatrixXd compute_step_matrix(const std::variant<double, std::vector<double>>& force_variant) const;

public:
    /**
     * @brief Sets stability parameter $\alpha \in (0, 2]$.
     * @param value Parameter value ($\alpha = 2.0$ gives standard Gaussian diffusion).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_alpha(double value);

    /**
     * @brief Sets skewness parameter $\beta \in [-1, 1]$.
     * @param value Parameter value (0 gives symmetric spread, negative leans left, positive leans right).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_beta(double value);

    /**
     * @brief Sets scale parameter $\sigma > 0$.
     * @param value Spread/dispersion factor of the diffusion profile.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_sigma(double value);

    /**
     * @brief Sets physical length of 1D spatial domain ($L$).
     * @param value Total domain length in simulation units.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_length(double value);

    /**
     * @brief Applies a constant, uniform external force across the domain.
     * @param value Constant force value.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_drift(double value);

    /**
     * @brief Applies custom force values provided at each grid point.
     * @param value List containing force values for all grid points.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_force(std::vector<double> value);

    /**
     * @brief Applies position-dependent force using a mathematical function $f(x)$.
     * @param value Function taking spatial position $x$ and returning local force.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_force(std::function<double(double)> value);

    /**
     * @brief Sets time step size ($\Delta t$).
     * @param value Duration of a single step forward in time.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_dt(double value);

    /**
     * @brief Sets the total number of grid intervals ($N$).
     * @param value Number of segments (creates $N + 1$ total grid points).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_num_intervals(Eigen::Index value);

    /**
     * @brief Sets time integration scheme parameter $\theta$.
     * @param value Value weighting scheme ($0.0$ explicit, $0.5$ semi-implicit, $1.0$ implicit).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_theta(double value);

    /**
     * @brief Sets detail level for terminal logging output.
     * @param value Verbosity level for built objects.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_verbosity(int value);

    /**
     * @brief Sets column delimiter string used when saving output data.
     * @param value Delimiter sequence (e.g., tab or comma).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_delimiter(std::string value);

    /**
     * @brief Sets completion progress percentage required between log updates.
     * @param value Progress update threshold percentage (%).
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_log_interval_percent(double value);

    /**
     * @brief Sets initial state values using an Eigen vector.
     * @param value Vector holding initial values for all grid points.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_initial_conditions(const Eigen::VectorXd &value);

    /**
     * @brief Sets initial state values using a standard vector.
     * @param value List holding initial values for all grid points.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_initial_conditions(const std::vector<double> &value);

    /**
     * @brief Places origin $x = 0.0$ at a specific grid point index.
     * @param value Grid index representing position $x = 0.0$.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_zero_index(Eigen::Index value);

    /**
     * @brief Places origin $x = 0.0$ at a physical distance from the left boundary.
     * @param value Distance value in length units from the left edge.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_zero_distance(double value);

    /**
     * @brief Places origin $x = 0.0$ at an exact center of domain.
     * @return SchemerBuilder& Reference to this builder for chaining.
     */
    SchemerBuilder &set_zero_middle();

    /**
     * @brief Builds and returns a fully configured simulation runner ready to execute.
     * @return SchemerRunner Ready-to-run simulation solver object.
     */
    [[nodiscard]] SchemerRunner build() const;

    /**
     * @brief Prepares and returns a complete simulation setup structure.
     * @return SimulationSetup Prepared setup with configuration, matrices, and initial states.
     */
    [[nodiscard]] SimulationSetup build_setup() const;
};