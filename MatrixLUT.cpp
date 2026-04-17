#include "MatrixLUT.h"

double MatrixLUT::lambda(Index k, Index m, std::vector<double> (&cacheRef)[2]) {
    if (k < 0) return 0.0;

    Index i = 0;
    if (k == i) return 1.0;

    if (k == 0 || k == this->size - 1) {
        if (std::isnan(cacheRef[0][m]))
            cacheRef[0][m] = pow(m-1, gamma+1) - (m - gamma - 1) * pow(m, gamma);
        return cacheRef[0][m];
    }

    if (std::isnan(cacheRef[1][m]))
        cacheRef[1][m] = pow(m+1, gamma+1) - 2*pow(m, gamma+1) + pow(m-1, gamma+1);

    return cacheRef[1][m];
}

MatrixXd MatrixLUT::calculateR() {
    MatrixXd R = MatrixXd::Zero(size, size);
    if (alpha == 1) {
        for (Index m = 0; m < size; m++)
            R(m, m) = 2; //lambda(0)

        for (Index n = 1; n < size-1; n++)
            for (Index m = 0; m < size - n; m++)
                R(m+n, m) = 2.0/(2.0*n + 1.0) - 2.0/(2.0*(n-1) + 1.0);

        return R;
    }

    std::vector <double> cache[2];
    for (auto & i : cache) {
        i.resize(size);
        i.assign(size, std::numeric_limits<double>::quiet_NaN());
    }

    if (alpha < 1) {
        for (Index m = 0; m < size; m++)
            R(m, m) = -3.0/2.0 * lambda(m, 0, cache);

        for (Index m = 0; m < size-1; m++)
            R(m+1, m) = -3.0/2.0 * lambda(m, 1, cache) + 2.0 * lambda(m+1, 0, cache);

        for (Index n = 2; n < size-1; n++)
            for (Index m = 0; m < size - n; m++)
                R(m+n, m) = -3.0/2.0*lambda(m, n, cache) + 2.0*lambda(m+1, n-1, cache) - 1.0/2.0*lambda(m+2, n-2, cache);
    }

    if (alpha > 1) {
        for (Index m = 0; m < size; m++)
            R(m, m) = 1.0 / 12.0 * lambda(m - 1, 2, cache) - 4.0 / 3.0 * lambda(m, 1, cache) + 5.0 / 2.0 * lambda(m + 1, 0, cache);

        for (Index m = 0; m < size-1; m++)
            R(m, m+1) = 1.0/12.0 * lambda(m, 1, cache) - 4.0/3.0 * lambda(m+1, 0, cache);

        for (Index m = 0; m < size-2; m++)
            R(m, m+1) = 1.0/12.0 * lambda(m+1, 0, cache);

        for (Index m = 0; m < size-1; m++)
            R(m, m + 1) = 1.0 / 12.0 * lambda(m - 1, 3, cache) - 4.0 / 3.0 * lambda(m, 2, cache) + 5.0 / 2.0 *
                          lambda(m + 1, 1, cache) - 4.0 / 3.0 * lambda(m + 2, 0, cache);

        for (Index n = 2; n < size-1; n++)
            for (Index m = 0; m < size - n; m++)
                R(m + n, m) = 1.0 / 12.0 * lambda(m - 1, n+2, cache) - 4.0 / 3.0 * lambda(m, n+1, cache) + 5.0 / 2.0 *
                              lambda(m + 1, n, cache) - 4.0 / 3.0 * lambda(m + 2, n-1, cache) + 1.0/12.0 * lambda(m + 3, n-2, cache);
    }

    return R;
}

MatrixXd MatrixLUT::calculateV() {
    MatrixXd result = MatrixXd::Zero(size, size);
    for (Index m = 0; m < size; ++m)
        V(m, m) = 3;

    for (Index m = 0; m < size - 1; ++m)
        V(m + 1, m) = -4;

    for (Index m = 0; m < size - 2; ++m)
        V(m + 2, m) = 1;
    return result;
}


MatrixLUT::MatrixLUT(const Index size, const double alpha) {
    this->size = size;
    this->alpha = alpha;
    this->gamma = std::ceil(alpha) - alpha;

    this->Id = MatrixXd::Identity(size, size);

    this->V = calculateV();

    MatrixXd R = calculateR();
    this->A = R + R.transpose();
}
