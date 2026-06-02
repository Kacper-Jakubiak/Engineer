#pragma once
#include <Eigen/Dense>
#include <variant>

struct Params {
    double alpha;
    double beta;
    double sigma;
    double length;
    std::variant<double, std::vector<double> > force_mode;

    double dt;
    Eigen::Index num_intervals;
    double theta;

    std::string delimiter;
    double log_interval_percent;
    int verbose;

    Eigen::VectorXd initial_values;
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

    void initialize_params();

    void initialize_matrices();

    [[nodiscard]] std::vector<double> get_lambdas() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;

    [[nodiscard]] Eigen::MatrixXd get_drift() const;

    [[nodiscard]] Eigen::MatrixXd get_force() const;

public:
    static constexpr double ALPHA_EPSILON = 1e-8;

    explicit Schemer(Params params);

    void reset();

    void run(int steps, int save_every = 0, std::ostream *os = nullptr);

    void save_result(const std::string &filepath) const;

    Schemer(const Schemer &) = delete;

    Schemer &operator=(const Schemer &) = delete;
};
