#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>

Schemer::Schemer(double alpha, double beta, double K, double dt, Eigen::Index I, double length, double mi,
                 double theta, int verbose, Eigen::VectorXd initial_values, Eigen::VectorXd force) : alpha(alpha), beta(beta), K(K), dt(dt), I(I), length(length),
                                              theta(theta), verbose(verbose), mi(mi), initial_values(std::move(initial_values)), force(std::move(force)) {
    initialize_params();
    initialize_matrices();

    reset_simulation();
}

void Schemer::reset_simulation() {
    history.clear();
    history.push_back(initial_values);
    current = initial_values;
}

void Schemer::run(const int steps) {
    const int log_interval = std::max(1, steps / 20);
    for (int h = 0; h < steps; h++) {
        current = step_matrix * current;
        history.push_back(current);
        if (verbose > 1 && h % log_interval == 0) {
            std::cout << "Progress: " << (100.0 * h / steps) << "%\n";
        }
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
