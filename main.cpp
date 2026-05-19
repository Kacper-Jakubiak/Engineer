#include <iostream>
#include <format>
#include <fstream>

#include "include/Schemer.h"
#include "include/SchemerBuilder.h"
using namespace std;

int main(const int argc, char *argv[]) {
    double alpha, length, dt, K, mi, theta, beta;
    int verbose, steps;
    Eigen::Index I;

    if (argc < 10) {
        cout << "Enter: I length dt K mi theta steps verbose alpha\n";
        cin >> I >> length >> dt >> K >> mi >> theta >> steps >> verbose >> beta >> alpha;
    } else {
        I = std::stoi(argv[1]);
        length = std::stod(argv[2]);
        dt = std::stod(argv[3]);
        K = std::stod(argv[4]);
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
        std::cout << "I " << I << "\nalpha " << alpha << "\nlength " << length << "\nbeta " << beta << "\ndt " << dt << "\nK " << K <<
                "\nmi " << mi << "\ntheta " << theta << "\nsteps " << steps << "\nverbose " << verbose << '\n';

    std::string filename;
    cout << "Enter filename: ";
    // cin >> filename;
    std::ofstream out_stream(std::string(PROJECT_ROOT) + "/history.txt");

    cout << "BUILDING..." << std::endl;
    Schemer simulator = SchemerBuilder()
            .set_alpha(alpha)
            .set_beta(beta)
            .set_sigma(K)
            .set_dt(dt)
            .set_grid_size(I)
            .set_length(length)
            .set_drift(mi)
            .set_theta(theta)
            .set_verbosity(verbose)
            .set_zero_point(I/2)
            // .set_mi(10.0)
            .set_force([](double x){return -x;})
            .build();
    cout << "RUNNING..." << std::endl;
    simulator.run(steps, 1, &out_stream);
    cout << "FINISHED." << std::endl;
    // simulator.save_result(filename + ".txt");
}
