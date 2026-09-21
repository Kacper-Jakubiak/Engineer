#include <iostream>
#include <format>
#include <fstream>

#include "include/Schemer.h"
#include "capala/hist.h"
using namespace std;

double (*force_from_int(const int n))(const double x) {
    switch (n) {
        case 0:
            return [](const double x) { return 0.0; };
        case 1:
            return [](const double x) { return -2 * x; };
        case 2:
            return [](const double x) { return -4.0 * std::sin(x); };
        case 3:
            return [](const double x) { return std::max(-50.0, std::min(x - x * x * x, 50.0)); };
        case 4:
            return [](const double x) { return -(4.0 * x) / (1.0 + x * x); };
        default:
            throw std::invalid_argument("Invalid force input");
    }
}


int main(const int argc, char *argv[]) {
    Config config;
    int force_input;
    int steps;


    if (argc >= 9) {
        config.alpha = std::stod(argv[1]);
        config.beta = std::stod(argv[2]);
        config.sigma = std::stod(argv[3]);
        config.length = std::stod(argv[4]);
        config.dt = std::stod(argv[5]);
        config.num_intervals = std::stoi(argv[6]);
        force_input = std::stoi(argv[7]);
        steps = std::stoi(argv[8]);
    } else {
        cout << "Enter: alpha, beta, sigma, length" << std::endl;
        cin >> config.alpha >> config.beta >> config.sigma >> config.length;
        cout << "Enter: dt, num_intervals" << std::endl;
        cin >> config.dt >> config.num_intervals;
        cout << "Enter: force (1, 2, 3, 4)" << std::endl;
        cin >> force_input;
        cout << "Enter: steps" << std::endl;
        cin >> steps;
    }

    // GŁÓWNY MAIN

    const std::string result_filepath = std::string(PROJECT_ROOT) + "/result.csv";
    const std::string history_filepath = std::string(PROJECT_ROOT) + "/history.csv";
    std::ofstream history_stream(history_filepath);

    auto force = force_from_int(force_input);

    cout << "BUILDING..." << std::endl;
    SchemerRunner simulator = SchemerBuilder()
            .set_alpha(config.alpha)
            .set_beta(config.beta)
            .set_sigma(config.sigma)
            .set_length(config.length)
            .set_dt(config.dt)
            .set_num_intervals(config.num_intervals)
            .set_log_interval_percent(1)
            .set_force(force)
            .build();


    cout << "RUNNING..." << std::endl;
    simulator.run(steps, &history_stream, 50);

    cout << "FINISHED." << std::endl;
    simulator.save_result(result_filepath);

    // TRAJEKTORIE

    return 0;

    cout << "TRAJECTORIES..." << std::endl;
    const int less_steps = steps / 10;
    const auto histogram = get_histogram(config.alpha, config.beta, config.sigma, config.length,
                                         static_cast<int>(config.num_intervals), (config.dt * steps) / less_steps,
                                         less_steps, force);

    std::ofstream histogram_stream(std::string(PROJECT_ROOT) + "/histogram.csv");
    for (const double v: histogram) {
        histogram_stream << v << "\t";
    }
    histogram_stream << std::endl;
}
