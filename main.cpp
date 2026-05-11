#include <iostream>
#include <format>

#include "Schemer.h"
using namespace std;

int main(int argc, char *argv[]) {
    double alpha, length, dt, K, mi, theta;
    int verbose;
    Eigen::Index I, J;

    if (argc < 9) {
        cout << "Enter: I length dt K mi theta J verbose alpha\n";
        cin >> I >> length >> dt >> K >> mi >> theta >> J >> verbose >> alpha;
    } else {
        I = std::stoi(argv[1]);
        length = std::stod(argv[2]);
        dt = std::stod(argv[3]);
        K = std::stod(argv[4]);
        mi = std::stod(argv[5]);
        theta = std::stod(argv[6]);
        J = std::stoi(argv[7]);
        verbose = std::stoi(argv[8]);
        if (argc == 9) {
            cout << "Enter alpha: " << endl;
            cin >> alpha;
        } else alpha = std::stod(argv[9]);
    }
    if (verbose > 0)
        std::cout << "I " << I << "\nalpha " << alpha << "\nlength " << length << "\ndt " << dt << "\nK " << K <<
                "\nmi " << mi << "\ntheta " << theta << "\nJ " << J << "\nverbose " << verbose << '\n';

    Schemer simulator{alpha, 0.0, K, dt, I, length, mi, theta, verbose};
    simulator.run(J);
    simulator.save_history("result_" + std::format("{:.2f}", alpha) + ".txt");
}
