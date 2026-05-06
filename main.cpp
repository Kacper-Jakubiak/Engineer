#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

using namespace Eigen;
using namespace std;
Index I;

void save_history(const vector<VectorXd> &history, const string &filename) {
    IOFormat fmt(
        4,
        DontAlignCols,
        ";",
        "\n",
        ""
    );
    const string path = std::string(PROJECT_ROOT) + "/" + filename;
    ofstream out(path);
    for (const VectorXd &v: history) {
        out << v.transpose().format(fmt) << endl;
    }
    cout << "Saved to " << path << endl;
}

void save_matrix(const MatrixXd &matrix, const string &filename) {
    ofstream out(std::string(PROJECT_ROOT) + "/" + filename);
    out << std::fixed << std::setprecision(3) << std::setw(5);
    out << matrix << endl;
    cout << "Saved " + filename << endl;
    out.close();
}

double get_lambda(Index k, Index m, double gamma) {
    if (k < 0 || m <= 0 || I < k)
        throw std::runtime_error("Invalid lambda invocation");

    // if (k == 0 || k == I) {
    //     return pow(m - 1, gamma + 1) - (m - gamma - 1) * pow(m, gamma);
    // }

    return pow(m + 1, gamma + 1) - 2 * pow(m, gamma + 1) + pow(m - 1, gamma + 1);
}

double get_lambda(Index m) {
    return 2.0 / (1.0 + static_cast<double>(m));
}

MatrixXd calculate_V(const double sign) {
    MatrixXd V = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        V(i, i) = 3.0;
        if (i >= 1)
            V(i, i - 1) = -4.0;
        if (i >= 2)
            V(i, i - 2) = 1.0;
    }
    if (sign < 0)
        return -V.transpose();
    return V;
}

MatrixXd calculate_R1(double gamma) {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, i - k, gamma);
            R(i, k) += 3.0 * lambda;
            if (k >= 1)
                R(i, k - 1) += -4.0 * lambda;
            if (k >= 2)
                R(i, k - 2) += 1.0 * lambda;
        }
    }
    return -1.0 / 2 * R;
}

MatrixXd calculate_R2(double gamma) {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, i - k, gamma);
            if (k + 2 <= I)
                R(i, k + 2) += -1.0 * lambda;
            if (k + 1 <= I)
                R(i, k + 1) += 16.0 * lambda;
            R(i, k) += -30.0 * lambda;
            if (k >= 1)
                R(i, k - 1) += 16.0 * lambda;
            if (k >= 2)
                R(i, k - 2) += -1.0 * lambda;
        }
    }
    return -1.0 / 12 * R;
}

MatrixXd calculate_R3() {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = 0; k <= i; k++) {
            const double lambda = get_lambda(i - k);
            R(i, k) += 1.0 * lambda;
            if (k >= 1)
                R(i, k - 1) += -1.0 * lambda;
        }
    }
    return -1.0 * R;
}

MatrixXd calculate_R3_ant() {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index n = 0; n < I + 1; n++) {
        const double lambda = get_lambda(n);
        for (Index m = 0; m + n < I + 1; m++) {
            R(m + n, m) += 1.0 * lambda;
            if (m + n + 1 < I + 1)
                R(m + n + 1, m) += -1.0 * lambda;
        }
    }
    return -1.0 * R;
}

MatrixXd get_R(int state, double gamma) {
    switch (state) {
        case 1:
            return calculate_R1(gamma);
        case 2:
            return calculate_R2(gamma);
        case 3:
            return calculate_R3();
        default:
            throw std::invalid_argument("Invalid state");
    }
}

MatrixXd calculate_R1t(double gamma) {
    MatrixXd R = MatrixXd::Zero(I + 1, I + 1);
    for (Index i = 0; i < I + 1; i++) {
        for (Index k = i; k <= I; k++) {
            const double lambda = (k == i) ? 1.0 : get_lambda(k, k - i, gamma);
            R(i, k) += -3.0 * lambda;
            if (k + 1 <= I)
                R(i, k + 1) += 4.0 * lambda;
            if (k + 2 <= I)
                R(i, k + 2) += -1.0 * lambda;
        }
    }
    //Chyba tak?
    return 1.0 / 2 * R;
}

int main(int argc, char *argv[]) {
    IOFormat fmt(
        2,
        DontAlignCols,
        ",",
        "\n",
        "[", "]"
    );

    cout << std::fixed;

    double alpha, length, dt, K, mi, theta;
    int J, verbose;

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

    int state;
    if (std::abs(alpha - 1.0) < 1e-6)
        state = 3;
    else if (alpha <= 1.0)
        state = 1;
    else
        state = 2;

    if (verbose > 0) cout << "state: " << state << endl;
    constexpr double pi = std::numbers::pi;
    double gamma = std::ceil(alpha) - alpha;
    double dx = length / static_cast<double>(I);

    double omega;
    if (state != 3) {
        const double denominator = 2.0 * std::cos(pi * alpha / 2.0) * std::tgamma(2.0 + gamma) *
                                   std::pow(dx, alpha);
        omega = -K * dt / denominator;
    } else {
        omega = -K * dt / (pi * dx);
    }
    double ni = mi * dt / (2.0 * dx);

    if (verbose > 0)
        cout << "gamma: " << gamma << "\nomega: " << omega << "\nni: " << ni << endl;

    VectorXd f = VectorXd::Zero(I + 1);
    f(I / 2) = length / dx;
    const MatrixXd Id = MatrixXd::Identity(I + 1, I + 1);
    const MatrixXd V = calculate_V(ni);
    const MatrixXd R = get_R(state, gamma);
    const MatrixXd Rt = R.transpose().eval(); //calculate_R1t(gamma);

    const MatrixXd A = R + Rt;
    const MatrixXd lhs = Id + theta * (ni * V + omega * A);
    const MatrixXd rhs = Id - (1.0 - theta) * (ni * V + omega * A);

    vector<VectorXd> history;
    history.push_back(f);
    LDLT<MatrixXd> solver(lhs);
    for (int h = 0; h < J; h++) {
        if (verbose > 1) {
            cout << f.transpose().format(fmt) << endl;
            cout << f.sum() << endl;
        }
        const VectorXd next = solver.solve(rhs * f);
        f = next;
        history.push_back(f);
    }
    cout << f.transpose().format(fmt) << endl;
    std::stringstream ss;
    ss << I << "_" << length << "_" << dt << "_" << K << "_" << mi << "_" << theta << "_" << J << "_" << alpha << ".txt";
    save_history(history, ss.str());
}
