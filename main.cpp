#include <iostream>
#include <format>

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
    cin >> filename;

    cout << "BUILDING.." << std::endl;
    Schemer simulator = SchemerBuilder()
            .set_alpha(alpha)
            .set_beta(beta)
            .set_K(K)
            .set_dt(dt)
            .set_I(I)
            .set_length(length)
            .set_mi(mi)
            .set_theta(theta)
            .set_verbose(verbose)
            .build();
    cout << "RUNNING.." << std::endl;
    simulator.run(steps);
    cout << "FINISHED.." << std::endl;
    simulator.save_history(filename + ".txt");
}
