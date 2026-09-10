/**
 * @file SchemerRunner.h
 * @brief Header for the Schemer solver and its settings.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include "RunnerSetup.h"

/**
 * @class SchemerRunner
 * @brief Main class that runs the simulation over time.
 */
class SchemerRunner {
private:
    /// File to save progress logs.
    static constexpr std::string_view LOG_FILE_PATH = "progress.log";

    RunnerSetup system;                 ///< Saved settings.
    int steps_taken = 0;                   ///< Number of time steps done.
    Eigen::VectorXd current_values;        ///< Current values on the grid.

    /**
     * @brief Prints progress to a stream.
     * @param os The stream to write to.
     * @param progress_percent How much is done (in percent).
     */
    static void log(std::ostream &os, double progress_percent);

    /**
     * @brief Saves settings to a stream.
     * @param os The stream to write to.
     */
    void save_parameters(std::ostream &os) const;

    /**
     * @brief Saves the current grid values to a stream.
     * @param os The stream to write to.
     */
    void save_current_values(std::ostream &os) const;

public:
    /**
     * @brief Creates a Schemer object with given settings and starting values.
     * @param system
     */
    explicit SchemerRunner(RunnerSetup system);

    ~SchemerRunner() = default;

    SchemerRunner(const SchemerRunner&) = default;
    SchemerRunner& operator=(const SchemerRunner&) = default;

    SchemerRunner(SchemerRunner&&) = default;
    SchemerRunner& operator=(SchemerRunner&&) = default;

    /**
     * @brief Resets the simulation back to the start (time zero).
     */
    void reset();

    /**
     * @brief Runs the simulation for a number of steps.
     * @param steps How many time steps to run.
     * @param history_stream Optional stream to save steps along the way.
     * @param save_every How often to save to the history stream (0 means don't save).
     */
    void run(int steps, std::ostream *history_stream = nullptr, int save_every = 0);

    /**
     * @brief Saves the final result to a file.
     * @param filepath The file to save to.
     */
    void save_result(const std::string &filepath) const;
};