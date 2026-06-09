#pragma once
#include <Eigen/Dense>
#include <variant>

struct Params {
    double alpha = 0.0;
    double beta = 0.0;
    double sigma = 1.0;
    double length = 40.0;
    std::variant<double, std::vector<double> > force_mode = 0.0;

    double dt = 0.001;
    Eigen::Index num_intervals = 1000;
    double theta = 0.5;

    std::string delimiter = "\t";
    double log_interval_percent = 5.0;
    int verbose = 1;
};

class Schemer {
private:
    static constexpr std::string_view LOG_FILE_PATH = "progress.log";

    enum class ForceType {
        Force,
        Drift
    };

    enum class SchemeType {
        Subdiffusive,
        Superdiffusive,
        Cauchy
    };

    const Params params;
    const Eigen::VectorXd initial_values;
    int steps_taken;

    Eigen::MatrixXd step_matrix;
    Eigen::VectorXd current_values;
    Eigen::Index size;
    double dx;
    double L;
    double R;
    double omega;
    double n;
    SchemeType state;
    ForceType force_type;

    static void log(std::ostream &os, double progress_percent);

    void save_parameters(std::ostream &os) const;

    void save_current_values(std::ostream &os) const;

    void initialize_values();

    void initialize_matrices();

    [[nodiscard]] std::vector<double> get_lambdas() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;

    [[nodiscard]] Eigen::MatrixXd get_drift() const;

    [[nodiscard]] Eigen::MatrixXd get_force() const;

public:
    static constexpr double ALPHA_EPSILON = 1e-8;

    Schemer(Params params, Eigen::VectorXd initial_values);

    void reset();

    void run(int steps, std::ostream *history_stream = nullptr, int save_every = 0);

    void save_result(const std::string &filepath) const;

    Schemer(const Schemer &) = delete;

    Schemer &operator=(const Schemer &) = delete;
};
