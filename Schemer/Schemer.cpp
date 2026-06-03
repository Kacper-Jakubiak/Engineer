#include "../include/Schemer.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <Eigen/Dense>
#include <utility>
#include <chrono>
#include <iomanip>

Schemer::Schemer(Params params): params(std::move(params)) {
    initialize_values();
    initialize_matrices();
    reset();
}

void Schemer::log(std::ostream &os, const double progress_percent) {
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);

    os << '['
            << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
            << "] Progress: "
            << progress_percent
            << "%" << std::endl;
}

void Schemer::save_parameters(std::ostream &os) const {
    os << '#' << "alpha: " << params.alpha << '\n';
    os << '#' << "beta: " << params.beta << '\n';
    os << '#' << "sigma: " << params.sigma << '\n';
    os << '#' << "length: " << params.length << '\n';

    os << '#' << "dt: " << params.dt << '\n';
    os << '#' << "num_intervals: " << params.num_intervals << '\n';
    os << '#' << "theta: " << params.theta << '\n';

    os << '#' << "delimiter: " << params.delimiter << '\n';
    os << '#' << "log_interval_percent: " << params.log_interval_percent << '\n';
    os << '#' << "verbose: " << params.verbose << '\n';

    os << '#' << "initial_values: ";
    for (Eigen::Index i = 0; i < params.initial_values.size(); i++)
        os << params.initial_values(i) << params.delimiter;
    os << '\n';

    switch (force_type) {
        case ForceType::Drift: {
            const double mi = std::get<double>(params.force_mode);
            os << '#' << "force_mode: " << mi << std::endl;
            break;
        }
        case ForceType::Force: {
            const auto &force_vector = std::get<std::vector<double>>(params.force_mode);
            os << '#' << "force_mode: ";
            for (const double force_value : force_vector)
                os << force_value << params.delimiter;
            os << std::endl;
            break;
        }
        default:
            throw std::runtime_error("Invalid Force Type");
    }

    // if (force_type == ForceType::Drift) {
    //     const double mi = std::get<double>(params.force_mode);
    //     os << '#' << "force_mode: " << mi << '\n';
    // }
    // else if (const std::vector<double>* force_mode_2 = std::get_if<std::vector<double> >(&params.force_mode)) {
    //     os << '#' << "force_mode: ";
    //     for (size_t i = 0; i < (*force_mode_2).size(); i++)
    //         os << (*force_mode_2)[i] << params.delimiter;
    //     os << '\n';
    // }
}

void Schemer::save_current_values(std::ostream &os) const {
    for (Eigen::Index i = 0; i < current_values.size(); i++)
        os << current_values(i) << params.delimiter;
    os << std::endl;
}

void Schemer::reset() {
    current_values = params.initial_values;
}

void Schemer::run(const int steps, std::ostream *history_stream, const int save_every) {
    if (steps < 0)
        throw std::invalid_argument("amount of steps cannot be negative");
    if (save_every < 0)
        throw std::invalid_argument("save_every cannot be negative");
    if (save_every > 0 && history_stream == nullptr)
        std::cerr << "[WARNING] history_stream not provided, history will not be saved";

    const bool should_log_history = history_stream != nullptr && save_every > 0;
    if (should_log_history) {
        if (history_stream->tellp()<= 0)
            save_parameters(*history_stream);
    }

    std::ofstream log_stream;
    bool should_log_progress = params.log_interval_percent > 0;
    if (should_log_progress) {
        log_stream.open(LOG_FILE_PATH.data());
        if (!log_stream) {
            std::cerr << "[WARNING] Failed opening log file " + std::string(LOG_FILE_PATH) << "logs will not be saved";
            should_log_progress = false;
        }
    }
    const int log_interval = std::max(1, static_cast<int>(steps * params.log_interval_percent / 100.0));

    for (int step = 0; step < steps; step++) {
        if (should_log_progress && step % log_interval == 0)
            log(log_stream, 100.0 * step / steps);

        if (should_log_history && step % save_every == 0)
            save_current_values(*history_stream);

        current_values = step_matrix * current_values;
    }

    if (should_log_history)
        save_current_values(*history_stream);

    if (should_log_progress)
        log(log_stream, 100.0);
}


void Schemer::save_result(const std::string &filepath) const {
    std::ofstream out_stream(filepath);
    if (!out_stream) {
        throw std::runtime_error("Failed to open save file: " + filepath);
    }
    save_parameters(out_stream);
    save_current_values(out_stream);
    if (params.verbose > 0)
        std::cout << "Saved to " << filepath << std::endl;
}