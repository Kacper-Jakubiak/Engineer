#pragma once
#include <Eigen/Dense>
#include <variant>

struct PhysicsParams {
    double alpha;
    double beta;
    double sigma;
    double length;
    std::variant<double, std::vector<double>> force_mode;
};

struct SolverParams {
    double dt;
    Eigen::Index grid_points;
    double theta;
};

class Schemer {
private:
    enum class ForceType {
        Force,
        Drift
    };

    enum class SchemeType {
        Subdiffusive,
        Superdiffusive,
        Cauchy
    };

    PhysicsParams physics;
    SolverParams solving;
    const int verbose;
    const Eigen::VectorXd initial_values;


    Eigen::MatrixXd step_matrix;
    Eigen::VectorXd current;
    Eigen::Index size;
    double dx;
    double L;
    double R;
    double omega;
    double n;
    double mi;
    std::vector<double> force;
    SchemeType state;
    ForceType force_type;

    void initialize_params();

    void initialize_matrices();

    [[nodiscard]] std::vector<double> get_lambdas() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;
    [[nodiscard]] Eigen::MatrixXd get_drift() const;
    [[nodiscard]] Eigen::MatrixXd get_force() const;

public:
    Schemer(PhysicsParams physics, SolverParams solving, Eigen::VectorXd initial_values, int verbose);

    void reset_simulation();
    void run(int steps, int save_every = 0, std::ostream* os = nullptr);
    void save_result(const std::string &filename) const;

    Schemer(const Schemer &) = delete;

    Schemer &operator=(const Schemer &) = delete;
};
