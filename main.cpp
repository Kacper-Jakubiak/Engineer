#include <iostream>
#include <format>
#include <fstream>

#include "include/Schemer.h"
#include "include/SchemerBuilder.h"
#include "capala/hist.h"
using namespace std;

double force(const double x) {
    return x - (x * x * x);
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

    std::string filename;
    // cout << "Enter filename: ";
    // cin >> filename;
    std::ofstream out_stream(std::string(PROJECT_ROOT) + "/history.txt");

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
            .set_force(force)
            .build();
    cout << "RUNNING..." << std::endl;
    simulator.run(steps, 10, &out_stream);
    cout << "FINISHED." << std::endl;
    simulator.save_result(std::string(PROJECT_ROOT) + "/result.txt");

    auto histogram = get_histogram();
    std::ofstream hist_stream(std::string(PROJECT_ROOT) + "/histogram.txt");
    for (double v : histogram) {
        hist_stream << v << ";";
    }
    hist_stream << std::endl;
}
