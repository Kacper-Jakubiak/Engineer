#include <iostream>
#include <format>

#include "Schemer.h"
using namespace std;

int main(int argc, char *argv[]) {
    double alpha, length, dt, K, mi, theta, beta;
    int verbose;
    Eigen::Index I, J;

    if (argc < 10) {
        cout << "Enter: I length dt K mi theta J verbose alpha\n";
        cin >> I >> length >> dt >> K >> mi >> theta >> J >> verbose >> beta >> alpha;
    } else {
        I = std::stoi(argv[1]);
        length = std::stod(argv[2]);
        dt = std::stod(argv[3]);
        K = std::stod(argv[4]);
        mi = std::stod(argv[5]);
        theta = std::stod(argv[6]);
        J = std::stoi(argv[7]);
        verbose = std::stoi(argv[8]);
        beta = std::stod(argv[9]);
        if (argc == 10) {
            cout << "Enter alpha: " << endl;
            cin >> alpha;
        } else alpha = std::stod(argv[10]);
    }
    if (verbose > 0)
        std::cout << "I " << I << "\nalpha " << alpha << "\nlength " << length << "\nbeta " << beta << "\ndt " << dt << "\nK " << K <<
                "\nmi " << mi << "\ntheta " << theta << "\nJ " << J << "\nverbose " << verbose << '\n';

    std::string filename;
    cout << "Enter filename: ";
    cin >> filename;

    Schemer simulator{alpha, beta, K, dt, I, length, mi, theta, verbose};
    simulator.run(J);
    simulator.save_history(filename + ".txt");
    // simulator.save_history("mibeta_" + std::format("{:.0f}", mi) + std::format("{:.2f}", beta)+ ".txt");
}
