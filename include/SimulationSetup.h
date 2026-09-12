/**
 * @file SimulationSetup.h
 * @brief Prepared simulation setup.
 */

#pragma once
#include "Config.h"
#include <Eigen/Dense>
#include <variant>
#include <vector>

/**
 * @struct SimulationSetup
 * @brief Complete setup ready to run.
 *
 * Contains all information needed to execute a simulation.
 */
struct SimulationSetup {
    /// Settings for the simulation
    Config config;

    /// Pre-computed operation matrix for time stepping
    Eigen::MatrixXd step_matrix;

    /// Initial values for all grid points
    Eigen::VectorXd initial_state;

    /// Starting position on the grid
    Eigen::Index starting_index;

    /// Force type: constant value or per-point values
    std::variant<double, std::vector<double>> force_variant;
};
