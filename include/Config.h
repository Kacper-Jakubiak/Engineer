/**
 * @file Config.h
 * @brief Configuration settings for simulation runs.
 */

#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

/**
 * @struct Config
 * @brief Holds all simulation settings, physical parameters, grid specs, and file output options.
 */
struct Config {
    /**
     * @brief Stability parameter ($\alpha$).
     *
     * Controls the order of the fractional derivative, restricted to $(0, 2]$.
     * A value of 2.0 gives standard Gaussian diffusion. Values less than 2.0
     * produce superdiffusion (Lévy flights) with heavy tails, making long jumps
     * more likely.
     */
    double alpha = 0.0;

    /**
     * @brief Skewness parameter ($\beta$).
     *
     * Controls the left-right bias of the diffusion, restricted to $[-1, 1]$.
     * A value of 0 gives a balanced, symmetric spread, negative values lean left,
     * and positive values lean right.
     */
    double beta = 0.0;

    /**
     * @brief Scale parameter ($\sigma$).
     *
     * Controls the overall spread or width of the distribution. Must be positive ($> 0$).
     */
    double sigma = 1.0;

    /**
     * @brief Size of the 1D spatial domain.
     *
     * Defines the length of the grid in simulation units.
     */
    double length = 40.0;

    /**
     * @brief Number of grid intervals ($N$).
     *
     * Divides the domain into $N$ segments, resulting in $N + 1$ total grid points.
     */
    Eigen::Index num_intervals = 1000;

    /**
     * @brief Time step size ($\Delta t$).
     *
     * The duration of a single step forward in time for numerical integration.
     */
    double dt = 0.001;

    /**
     * @brief Time stepping method ($\theta$).
     *
     * - $0.0$: Explicit scheme (Forward Euler)
     * - $0.5$: Semi-implicit scheme (Crank-Nicolson)
     * - $1.0$: Implicit scheme (Backward Euler)
     */
    double theta = 0.5;

    /**
     * @brief Controls how much detail is printed to the terminal.
     *
     * - $0$: Silent mode (no messages)
     * - $1$: Status updates
     */
    int verbose = 1;

    /**
     * @brief Log progress percentage.
     *
     * Completion percentage step between logging progress.
     */
    double log_interval_percent = 5.0;

    /**
     * @brief Column separator for output files.
     *
     * The string sequence used between data columns when saving results.
     */
    std::string delimiter = "\t";

    /**
     * @brief Calculates uniform grid step size ($\Delta x$).
     * @return double Grid step size.
     */
    [[nodiscard]] constexpr double get_dx() const noexcept {
        return length / static_cast<double>(num_intervals);
    }

    /**
     * @brief Gets the total number of grid points ($N + 1$).
     * @return Eigen::Index Total grid point count.
     */
    [[nodiscard]] constexpr Eigen::Index get_size() const noexcept {
        return num_intervals + 1;
    }

    /**
     * @brief Gets position coordinates for all grid points.
     * @param starting_index Grid point index corresponding to $x = 0.0$.
     * @return std::vector<double> List of spatial positions.
     */
    [[nodiscard]] std::vector<double> get_coordinates(const Eigen::Index starting_index = 0) const {
        const double dx = get_dx();
        const auto total_points = get_size();
        std::vector<double> coords(static_cast<std::size_t>(total_points));

        for (Eigen::Index i = 0; i < total_points; ++i) {
            coords[static_cast<std::size_t>(i)] = dx * static_cast<double>(i - starting_index);
        }
        return coords;
    }
};