#pragma once
#include <Eigen/Dense>

using namespace Eigen;

class MatrixLUT {
public:
    MatrixXd A;
    Index size;
    double alpha;
    double gamma;

    explicit MatrixLUT(Index size, double alpha);

    static MatrixXd calculate_V(Index size, double ni);

private:
    double lambda(Index k, Index m, std::vector<double> (&cacheRef)[2]);

    MatrixXd calculateR();
};
