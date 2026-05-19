#pragma once
#include <Eigen/Dense>

class Schemer {
private:
    enum class SchemeType {
        Subdiffusive,
        Superdiffusive,
        Cauchy
    };

    const double alpha;
    const double beta;
    const double K;
    const double dt;
    const Eigen::Index I;
    const double length;
    const double theta;
    const int verbose;
    const double mi;
    Eigen::MatrixXd step_matrix;
    const Eigen::VectorXd initial_values;
    const std::vector<double> force;
    Eigen::VectorXd current;

    double dx;
    double L;
    double R;
    double omega;
    double n;
    double ni;
    SchemeType state;

    void initialize_params();

    void initialize_matrices();

    [[nodiscard]] std::vector<double> get_lambdas() const;

    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;
    [[nodiscard]] Eigen::MatrixXd get_drift() const;
    [[nodiscard]] Eigen::MatrixXd get_force() const;

public:
    Schemer(double alpha, double beta, double K,
            double dt, Eigen::Index I, double length, double mi, double theta, int verbose, Eigen::VectorXd initial_values, std::vector<double> force);

    void reset_simulation();
    void run(int steps, int save_every = 0, std::ostream* os = nullptr);
    void save_result(const std::string &filename) const;

    Schemer(const Schemer &) = delete;

    Schemer &operator=(const Schemer &) = delete;
};
