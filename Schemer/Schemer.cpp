#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>
#include <utility>
#include <chrono>
#include <iomanip>

Schemer::Schemer(PhysicsParams physics, SolverParams solving, FrontParams front, Eigen::VectorXd initial_values)
    : physics(std::move(physics)), solving(std::move(solving)), front(std::move(front)), initial_values(std::move(initial_values)) {

    initialize_params();
    initialize_matrices();
    reset();
}

void Schemer::log(std::ostream* os, const double progress_percent)
{
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);

    *os << '['
        << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
        << "] Progress: "
        << progress_percent
        << "%\n";
}

void Schemer::save_current_values(std::ostream* os) const {
    for (Eigen::Index i = 0; i < current_values.size(); i++)
        *os << current_values(i) << front.delimiter;
    *os << std::endl;
}


void Schemer::reset() {
    current_values = initial_values;
}

void Schemer::run(const int steps, const int save_every, std::ostream* os) {
    if (steps <= 0)
        throw std::invalid_argument("steps must be > 0");
    if (save_every < 0)
        throw std::invalid_argument("save_every must be >= 0");
    if (save_every > 0 && os == nullptr)
        throw std::invalid_argument("stream is required");

    std::ofstream log_file(LOG_FILE_PATH.data());
    const int log_interval = std::max(1,  static_cast<int>(steps * front.log_interval_percent / 100.0));

    for (int step = 0; step < steps; step++) {
        if (front.log_interval_percent > 0.0 && step % log_interval == 0)
            log(&log_file, 100.0 * step / steps);

        if (save_every > 0 && step % save_every == 0)
            save_current_values(os);

        current_values = step_matrix * current_values;
    }

    if (save_every > 0)
        save_current_values(os);

    if (front.log_interval_percent > 0.0)
        log(&log_file, 100.0);
}


void Schemer::save_result(const std::string &filepath) const {
    std::ofstream out_stream(filepath);
    for (Eigen::Index i = 0; i < current_values.size(); i++)
        out_stream << current_values(i) << ";";
    out_stream << std::endl;
    if (front.verbose > 0)
        std::cout << "Saved to " << filepath << std::endl;
}
