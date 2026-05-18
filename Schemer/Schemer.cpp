#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>

Schemer::Schemer(double alpha, double beta, double K, double dt, Eigen::Index I, double length, double mi,
                 double theta, int verbose) : alpha(alpha), beta(beta), K(K), dt(dt), I(I), length(length),
                                              theta(theta), verbose(verbose), mi(mi) {
    initialize_params();
    initialize_matrices();

    initial = Eigen::VectorXd::Zero(I + 1);
    initial(I / 2) = length / dx;

    const Eigen::PartialPivLU<Eigen::MatrixXd> solver(Lhs);
    step_matrix = solver.solve(Rhs);

    reset_simulation();
}

void Schemer::reset_simulation() {
    history.clear();
    history.push_back(initial);
    current = initial;
}

void Schemer::run(const int steps) {
    for (int h = 0; h < steps; h++) {
        current = step_matrix * current;
        history.push_back(current);
    }
}


void Schemer::save_history(const std::string &filename) const {
    Eigen::IOFormat fmt(4, Eigen::DontAlignCols, ";", "\n", "");
    const std::string path = std::string(PROJECT_ROOT) + "/" + filename;
    std::ofstream out(path);
    for (const Eigen::VectorXd &v: history) {
        out << v.transpose().format(fmt) << std::endl;
    }
    std::cout << "Saved to " << path << std::endl;
}
