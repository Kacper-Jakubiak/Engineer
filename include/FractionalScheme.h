/**
 * @file FractionalScheme.h
 * @brief Builds operation matrices for simulation.
 */

#pragma once

#include "Config.h"
#include <Eigen/Dense>
#include <vector>
#include <variant>

/**
 * @class FractionalScheme
 * @brief Creates matrices for time stepping.
 *
 * Builds the matrices needed to advance the simulation through time.
 */
class FractionalScheme {
public:
    /// Tolerance for detecting alpha = 0
    static constexpr double ALPHA_EPSILON = 1e-8;

    /**
     * @brief Initialize with settings.
     * @param config The simulation settings.
     */
    explicit FractionalScheme(const Config& config);
    
    /**
     * @brief Build the main diffusion matrix.
     * @return Eigen::MatrixXd The diffusion matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_diffusion_matrix() const;

    /**
     * @brief Build the force matrix.
     * @param force_mode Either constant force or force per grid point.
     * @return Eigen::MatrixXd The force matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_force_matrix(const std::variant<double, std::vector<double>>& force_mode) const;

private:
    const Config& config;

    double L;
    double R;
    double omega;
    double cauchyOmega;
    double n;

    /**
     * @brief Get weight coefficients.
     * @return std::vector<double> Weight values.
     */
    [[nodiscard]] std::vector<double> get_lambdas() const;

    /**
     * @brief Build matrix for case alpha < 1.
     * @return Eigen::MatrixXd The matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_M1() const;

    /**
     * @brief Build matrix for case alpha > 1.
     * @return Eigen::MatrixXd The matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_M2() const;

    /**
     * @brief Build matrix for case alpha = 1.
     * @return Eigen::MatrixXd The matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_M3() const;

    /**
     * @brief Build constant force matrix.
     * @param mi The force value.
     * @return Eigen::MatrixXd The force matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_drift(double mi) const;

    /**
     * @brief Build variable force matrix.
     * @param force Force values for each grid point.
     * @return Eigen::MatrixXd The force matrix.
     */
    [[nodiscard]] Eigen::MatrixXd build_force(const std::vector<double>& force) const;
};