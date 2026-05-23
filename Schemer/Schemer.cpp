#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>
#include <utility>

Schemer::Schemer(PhysicsParams physics, SolverParams solving, Eigen::VectorXd initial_values, const int verbose)
    : physics(std::move(physics)), solving(std::move(solving)), verbose(verbose), initial_values(std::move(initial_values)) {

    initialize_params();
    initialize_matrices();
    reset_simulation();
}


void Schemer::reset_simulation() {
    current = initial_values;
}

void Schemer::run(const int steps, const int save_every, std::ostream* os) {
    if (steps <= 0)
        throw std::invalid_argument("steps must be > 0");
    if (save_every < 0)
        throw std::invalid_argument("save_every must be >= 0");
    if (save_every > 0 && os == nullptr)
        throw std::invalid_argument("stream is required");

    const int log_interval = std::max(1, steps / 20);

    for (int h = 0; h < steps; h++) {
        if (verbose > 0 && h % log_interval == 0) {
            std::cout << "\rProgress: " << (100.0 * h / steps) << "%" << std::flush;
        }
        if (save_every > 0 && h % save_every == 0) {
            for (Eigen::Index i = 0; i < current.size(); i++)
                *os << current(i) << ";";
            *os << std::endl;
        }
        current = step_matrix * current;
    }
    if (save_every > 0) {
        for (Eigen::Index i = 0; i < current.size(); i++)
            *os << current(i) << ";";
        *os << std::endl;
    }
    if (verbose > 0) {
        std::cout << "\rProgress: " << 100.0 << "%" << std::endl;
    }
}


void Schemer::save_result(const std::string &filename) const {
    const std::string path = std::string(PROJECT_ROOT) + "/" + filename;
    std::ofstream out_stream(path);
    for (Eigen::Index i = 0; i < current.size(); i++)
        out_stream << current(i) << ";";
    out_stream << std::endl;
    if (verbose > 0)
        std::cout << "Saved to " << path << std::endl;
}
