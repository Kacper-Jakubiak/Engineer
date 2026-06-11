#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "alfaStabilny.h"
#include "hist.h"

#include <functional>
#define NUM_PARTICLES 100000

std::vector<double> generate_histogram(const std::vector<double> &data, int num_bins, double min_range,
                                       double max_range) {
    std::vector<int> counts(num_bins, 0);

    double length = max_range - min_range;
    double bin_width = length / num_bins;
    int out_of_range = 0;

    for (double value: data) {
        if (value < min_range || value >= max_range) {
            out_of_range++;
            continue;
        }

        const int bin_idx = static_cast<int>((value - min_range) / bin_width);

        if (bin_idx < 0 || bin_idx >= num_bins) {
            out_of_range++;
            continue;
        }

        counts[bin_idx]++;
    }

    std::cout << out_of_range << " / " << data.size() << " out of range";

    std::vector density(num_bins, 0.0);

    const double normalization_factor = static_cast<double>(data.size()) * bin_width;
    for (int i = 0; i < num_bins; ++i)
        density[i] = counts[i] / normalization_factor;

    return density;
}


std::vector<double> generate_positions(alfaStabilny &levyGenerator, int num_particles, double dt, int steps,
                                       double noise, const std::function<double(double)> &force) {
    vector<double> positions;
    positions.reserve(num_particles);

    int update_interval = std::max(1, num_particles / 100);

    for (int id = 0; id < num_particles; ++id) {
        if (id % update_interval == 0) {
            int progress = (id * 100) / num_particles;
            std::cout << "\rProgress: " << progress << "%" << std::flush;
        }

        double x = 0.0;

        for (int step = 0; step < steps; ++step) {
            double drift = force(x) * dt;

            double jump = levyGenerator.losuj() * noise;

            x += drift + jump;
        }
        positions.push_back(x);
    }
    cout << "\rProgress: 100%" << std::endl;

    return positions;
}


std::vector<double> get_histogram(double alpha, double beta, double sigma, double length, int num_intervals, double dt,
                                  int steps, const std::function<double(double)> &force) {
    std::random_device rd;
    std::mt19937 gen(rd());

    alfaStabilny levyGenerator(alpha, beta, 0.0, sigma, &gen);

    double noise_scaling = std::pow(dt, 1.0 / alpha);

    auto positions = generate_positions(levyGenerator, NUM_PARTICLES, dt, steps, noise_scaling, force);
    auto histogram = generate_histogram(positions, num_intervals + 1, -length / 2, length / 2);

    return histogram;
}
