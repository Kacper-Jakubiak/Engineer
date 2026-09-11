/**
 * @file Config.h
 * @brief Settings structure for simulations.
 */

#pragma once

#include <vector>
#include <Eigen/Dense>

/**
 * @struct Config
 * @brief Holds all simulation settings.
 */
struct Config {
    double alpha = 0.0;  ///< Main parameter
    double beta = 0.0;   ///< Secondary parameter
    double sigma = 1.0;  ///< Scale factor

    double length = 40.0;              ///< Domain size
    Eigen::Index num_intervals = 1000; ///< Number of grid points

    double dt = 0.001;  ///< Time step size
    double theta = 0.5; ///< Time stepping parameter

    int verbose = 1;                   ///< Output detail level
    double log_interval_percent = 5.0; ///< How often to log progress (%)
    std::string delimiter = "\t";      ///< Column separator in output

    /**
     * @brief Calculate grid spacing.
     * @return double Grid step size.
     */
    [[nodiscard]] double get_dx() const {
        return length / static_cast<double>(num_intervals);
    }

    /**
     * @brief Get the total number of grid points.
     * @return Eigen::Index Grid size.
     */
    [[nodiscard]] Eigen::Index get_size() const {
        return num_intervals + 1;
    }

    /**
     * @brief Get position values for all grid points.
     * @param starting_index Starting position on the grid.
     * @return std::vector<double> Grid positions.
     */
    [[nodiscard]] std::vector<double> get_coordinates(const Eigen::Index starting_index) const {
        const double dx = get_dx();
        std::vector<double> coords(get_size());
        for (Eigen::Index i = 0; i < get_size(); ++i) {
            coords[i] = dx * static_cast<double>(i - starting_index);
        }
        return coords;
    }
};