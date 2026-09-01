/**
 * @file Schemer.h
 * @brief Header for the Schemer solver and its settings.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include <variant>
#include <vector>

/**
 * @struct Params
 * @brief Settings for the solver.
 */
struct Params {
    double alpha = 0.0;        ///< Fractional order (stability).
    double beta = 0.0;         ///< Skewness (asymmetry).
    double sigma = 1.0;        ///< Diffusion rate (scale).
    double length = 40.0;      ///< Total size of the space.

    /// External force: a single constant speed (double) or a custom array of forces (vector).
    std::variant<double, std::vector<double>> force_mode = 0.0;

    Eigen::Index starting_index = 0; ///< The grid index where position is zero.

    double dt = 0.001;                 ///< Size of each time step.
    Eigen::Index num_intervals = 1000; ///< Number of grid pieces.
    double theta = 0.5;                ///< Weight for time steps.

    std::string delimiter = "\t";      ///< Separator for text files.
    double log_interval_percent = 5.0; ///< How often to print progress (in percent).
    int verbose = 1;                   ///< How much info to print to the screen.
};

/**
 * @class Schemer
 * @brief Main class that runs the simulation over time.
 */
class Schemer {
private:
    /// File to save progress logs.
    static constexpr std::string_view LOG_FILE_PATH = "progress.log";

    /**
     * @enum ForceType
     * @brief Type of force: custom array (force) or constant drift.
     */
    enum class ForceType {
        Force, ///< Custom force array for each point.
        Drift  ///< A single constant drift speed.
    };

    /**
     * @enum SchemeType
     * @brief Math method used, chosen based on the alpha value.
     */
    enum class SchemeType {
        Subdiffusive,   ///< For alpha < 1.
        Superdiffusive, ///< For alpha > 1.
        Cauchy          ///< For alpha = 1.
    };

    const Params params;                   ///< Saved settings.
    const Eigen::VectorXd initial_values;  ///< Starting grid values.
    int steps_taken = 0;                   ///< Number of time steps done.

    Eigen::MatrixXd step_matrix;   ///< Matrix used to move forward in time.
    Eigen::VectorXd current_values;///< Current values on the grid.
    Eigen::Index size;             ///< Total number of grid points.
    double dx;                     ///< Space between grid points.
    double L;                      ///< Internal math parameter L.
    double R;                      ///< Internal math parameter R.
    double omega;                  ///< Internal scale factor.
    double n;                      ///< Internal power parameter.
    std::vector<double> coordinates;///< Physical positions of the grid points.
    SchemeType state;              ///< The active math method.
    ForceType force_type;          ///< The active force type.

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

    /**
     * @brief Sets up the starting positions and values.
     */
    void initialize_values();

    /**
     * @brief Sets up the math matrices.
     */
    void initialize_matrices();

    /**
     * @brief Calculates math weights (lambdas).
     * @return std::vector<double> The calculated weights.
     */
    [[nodiscard]] std::vector<double> get_lambdas() const;

    /**
     * @brief Builds the diffusion matrix for alpha < 1.
     * @return Eigen::MatrixXd The calculated matrix.
     */
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;

    /**
     * @brief Builds the diffusion matrix for alpha > 1.
     * @return Eigen::MatrixXd The calculated matrix.
     */
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;

    /**
     * @brief Builds the diffusion matrix for alpha = 1.
     * @return Eigen::MatrixXd The calculated matrix.
     */
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;

    /**
     * @brief Builds the drift matrix.
     * @return Eigen::MatrixXd The drift matrix.
     */
    [[nodiscard]] Eigen::MatrixXd get_drift() const;

    /**
     * @brief Builds the custom force matrix.
     * @return Eigen::MatrixXd The force matrix.
     */
    [[nodiscard]] Eigen::MatrixXd get_force() const;

public:
    /// Tolerance for checking if alpha is equal to 1.0.
    static constexpr double ALPHA_EPSILON = 1e-8;

    /**
     * @brief Creates a Schemer object with given settings and starting values.
     * @param params Settings to use.
     * @param initial_values Starting grid values.
     */
    Schemer(Params params, Eigen::VectorXd initial_values);

    // Prevent copying to avoid expensive or unsafe state duplication.
    Schemer(const Schemer &) = delete;
    Schemer &operator=(const Schemer &) = delete;

    // Prevent moving.
    Schemer(Schemer &&) noexcept = delete;
    Schemer &operator=(Schemer &&) noexcept = delete;

    /// Default destructor.
    ~Schemer() = default;

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