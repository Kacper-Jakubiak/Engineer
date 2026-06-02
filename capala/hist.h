#pragma once
#include <functional>

std::vector<double> get_histogram(double alpha, double beta, double sigma, double length, int num_intervals, double dt, int steps, const std::function<double(double)>& force);
