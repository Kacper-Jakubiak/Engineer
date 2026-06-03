#include <iostream>
#include <format>
#include <fstream>

#include "include/Schemer.h"
#include "include/SchemerBuilder.h"
#include "capala/hist.h"
using namespace std;

double get_force(const double x) {
    return -4.0 * std::sin(x) + x / 100.0;
}


int main(const int argc, char *argv[]) {
    double alpha, length, dt, sigma, mi, theta, beta;
    int verbose, steps;
    Eigen::Index num_intervals;

    if (argc < 10) {
        cout << "Enter: num_intervals length dt sigma mi theta steps verbose alpha\n";
        cin >> num_intervals >> length >> dt >> sigma >> mi >> theta >> steps >> verbose >> beta >> alpha;
    } else {
        num_intervals = std::stoi(argv[1]);
        length = std::stod(argv[2]);
        dt = std::stod(argv[3]);
        sigma = std::stod(argv[4]);
        mi = std::stod(argv[5]);
        theta = std::stod(argv[6]);
        steps = std::stoi(argv[7]);
        verbose = std::stoi(argv[8]);
        beta = std::stod(argv[9]);
        if (argc == 10) {
            cout << "Enter alpha: " << endl;
            cin >> alpha;
        } else alpha = std::stod(argv[10]);
    }
    if (verbose > 0)
        std::cout << "num_intervals " << num_intervals << "\nalpha " << alpha << "\nlength " << length << "\nbeta " << beta << "\ndt " << dt <<
                "\nsigma " << sigma <<
                "\nmi " << mi << "\ntheta " << theta << "\nsteps " << steps << "\nverbose " << verbose << '\n';

    // GŁÓWNY MAIN

    const std::string result_filepath = std::string(PROJECT_ROOT) + "/result.txt";
    const std::string history_filepath = std::string(PROJECT_ROOT) + "/history.txt";
    std::ofstream history_stream(history_filepath);

    cout << "BUILDING..." << std::endl;
    Schemer simulator = SchemerBuilder()
            .set_alpha(alpha)
            .set_beta(beta)
            .set_sigma(sigma)
            .set_dt(dt)
            .set_num_intervals(num_intervals)
            .set_length(length)
            .set_drift(mi)
            .set_theta(theta)
            .set_verbosity(verbose)
            .set_force(get_force)
            .set_log_interval_percent(1)
            .build();
    cout << "RUNNING..." << std::endl;
    simulator.run(steps, &history_stream, 5);
    cout << "FINISHED." << std::endl;
    simulator.save_result(result_filepath);

    // TRAJEKTORIE

    const int less_steps = steps / 20;
    const auto histogram = get_histogram(alpha, beta, sigma, length,
        static_cast<int>(num_intervals), (dt * steps) / less_steps, less_steps, get_force);

    std::ofstream histogram_stream(std::string(PROJECT_ROOT) + "/histogram.txt");
    for (const double v : histogram) {
        histogram_stream << v << "\t";
    }
    histogram_stream << std::endl;
}
