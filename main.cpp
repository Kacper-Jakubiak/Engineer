#include <iostream>
#include <format>
#include <fstream>

#include "include/Schemer.h"
#include "capala/hist.h"
using namespace std;

// double get_force(const double x) {
//     return 0.0;
// }

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
    Params params;
    int force_input;
    int steps;


    if (argc >= 9) {
        params.alpha = std::stod(argv[1]);
        params.beta = std::stod(argv[2]);
        params.sigma = std::stod(argv[3]);
        params.length = std::stod(argv[4]);
        params.dt = std::stod(argv[5]);
        params.num_intervals = std::stoi(argv[6]);
        force_input = std::stoi(argv[7]);
        steps = std::stoi(argv[8]);
    } else {
        cout << "Enter: alpha, beta, sigma, length" << std::endl;
        cin >> params.alpha >> params.beta >> params.sigma >> params.length;
        cout << "Enter: dt, num_intervals" << std::endl;
        cin >> params.dt >> params.num_intervals;
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
            .set_alpha(params.alpha)
            .set_beta(params.beta)
            .set_sigma(params.sigma)
            .set_length(params.length)
            .set_dt(params.dt)
            .set_num_intervals(params.num_intervals)
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
    const auto histogram = get_histogram(params.alpha, params.beta, params.sigma, params.length,
                                         static_cast<int>(params.num_intervals), (params.dt * steps) / less_steps,
                                         less_steps, force);

    std::ofstream histogram_stream(std::string(PROJECT_ROOT) + "/histogram.csv");
    for (const double v: histogram) {
        histogram_stream << v << "\t";
    }
    histogram_stream << std::endl;
}
