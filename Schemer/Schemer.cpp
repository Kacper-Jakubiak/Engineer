#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>

Schemer::Schemer(double alpha, double beta, double K, double dt, Eigen::Index I, double length, double mi,
                 double theta, int verbose, Eigen::VectorXd initial_values, std::vector<double> force) : alpha(alpha), beta(beta), K(K), dt(dt), I(I), length(length),
                                              theta(theta), verbose(verbose), mi(mi), initial_values(std::move(initial_values)), force(std::move(force)) {
    initialize_params();
    initialize_matrices();

    reset_simulation();
}

void Schemer::reset_simulation() {
    current = initial_values;
}

void Schemer::run(const int steps, const int save_every) {
    if (steps <= 0)
        throw std::invalid_argument("steps must be > 0");
    if (save_every < 0)
        throw std::invalid_argument("save_every must be >= 0");

    const int log_interval = std::max(1, steps / 20);
    std::ofstream out_stream;
    if (save_every > 0) {
        out_stream.open(std::string(PROJECT_ROOT) + "/" + "history.txt");
    }

    for (int h = 0; h < steps; h++) {
        current = step_matrix * current;
        if (verbose > 0 && h % log_interval == 0) {
            std::cout << "Progress: " << (100.0 * h / steps) << "%\n";
        }
        if (save_every > 0 && h % save_every == 0) {
            for (Eigen::Index i = 0; i < current.size(); i++)
                out_stream << current(i) << ";";
            out_stream << std::endl;
        }
    }
}


void Schemer::save_result(const std::string &filename) const {
    const std::string path = std::string(PROJECT_ROOT) + "/" + filename;
    std::ofstream out_stream(path);
    for (Eigen::Index i = 0; i < current.size(); i++)
        out_stream << current(i) << ";";
    out_stream << std::endl;
    // std::cout << "Saved to " << path << std::endl;
}
