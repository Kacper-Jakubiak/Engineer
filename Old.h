#pragma once
#include <Eigen/Dense>

class OldSchemer {
private:
    double alpha;
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
    double gamma;
    int state;
    double ni;
    std::vector<Eigen::VectorXd> history;

    [[nodiscard]] double lambda(Eigen::Index n) const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M1() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M2() const;
    [[nodiscard]] Eigen::MatrixXd get_diffusion_M3() const;


    [[nodiscard]]double get_lambda(Eigen::Index k, Eigen::Index m);
    [[nodiscard]]double get_lambda(Eigen::Index m);
    [[nodiscard]] Eigen::MatrixXd get_drift() const;
    [[nodiscard]] Eigen::MatrixXd get_R();
    [[nodiscard]] Eigen::MatrixXd get_R2();
    [[nodiscard]]Eigen::MatrixXd get_R3();


public:
    OldSchemer(double alpha, double K,
            double dt, Eigen::Index I, double length, double mi, double theta, int verbose);
    void run(Eigen::Index J);
    void save_history(const std::string &filename) const;
};