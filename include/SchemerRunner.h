/**
 * @file SchemerRunner.h
 * @brief Runs simulations over time.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include "SimulationSetup.h"

/**
 * @class SchemerRunner
 * @brief Executes simulations.
 *
 * Manages simulation state and runs time steps.
 */
class SchemerRunner {
private:
    /// File for progress logs
    static constexpr std::string_view LOG_FILE_PATH = "progress.log";

    SimulationSetup system;
    int steps_taken = 0;
    Eigen::VectorXd current_values;

    /**
     * @brief Log progress to a stream.
     * @param os Output stream.
     * @param progress_percent Completion percentage.
     */
    static void log(std::ostream &os, double progress_percent);

    /**
     * @brief Save configuration to a stream.
     * @param os Output stream.
     */
    void save_parameters(std::ostream &os) const;

    /**
     * @brief Save current values to a stream.
     * @param os Output stream.
     */
    void save_current_values(std::ostream &os) const;

public:
    /**
     * @brief Create a simulation runner.
     * @param system Setup with configuration and matrices.
     */
    explicit SchemerRunner(SimulationSetup system);

    ~SchemerRunner() = default;

    SchemerRunner(const SchemerRunner&) = default;
    SchemerRunner& operator=(const SchemerRunner&) = default;

    SchemerRunner(SchemerRunner&&) = default;
    SchemerRunner& operator=(SchemerRunner&&) = default;

    /**
     * @brief Reset to initial state.
     */
    void reset();

    /**
     * @brief Run simulation for a number of time steps.
     *
     * @param steps Number of steps to run.
     * @param history_stream Optional stream to save history (nullptr to skip).
     * @param save_every Save interval (0 = never).
     */
    void run(int steps, std::ostream *history_stream = nullptr, int save_every = 0);

    /**
     * @brief Save final result to file.
     *
     * @param filepath Path to save file.
     */
    void save_result(const std::string &filepath) const;
};