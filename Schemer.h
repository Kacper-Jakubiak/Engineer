#pragma once
#include <Eigen/Dense>

class Schemer {
private:
    double alpha;
    double beta;
    double K;
    double dt;
    Eigen::Index I;
    double length;
    double theta;
    int verbose;
    double mi;

    double dx;
    double L;
    double R;
    double omega;
    double n;
    int state;
    double ni;
    std::vector<Eigen::VectorXd> history;

    [[nodiscard]] double lambda(Eigen::Index n) const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;
    [[nodiscard]] Eigen::MatrixXd get_drift() const;

public:
    Schemer(double alpha, double beta, double K,
            double dt, Eigen::Index I, double length, double mi, double theta, int verbose);
    void run(Eigen::Index J);
    void save_history(const std::string &filename) const;
};
