#pragma once
#include <Eigen/Dense>

using namespace Eigen;

class MatrixLUT {
public:
    MatrixXd V;
    MatrixXd Id;
    MatrixXd A;
    Index size;
    double alpha;
    double gamma;

    explicit MatrixLUT(Index size, double alpha);
private:
    double lambda(Index k, Index m, std::vector<double> (&cacheRef)[2]);
    MatrixXd calculateR();
    MatrixXd calculateV();
};
