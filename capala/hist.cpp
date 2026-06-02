#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "alfaStabilny.h"
#include "hist.h"

double get_force(double x) {
    return x - (x * x * x);
}

std::vector<double> generate_histogram(const std::vector<double>& data, int num_bins, double min_range, double max_range) {
    std::vector<int> counts(num_bins, 0);

    double bin_width = (max_range - min_range) / num_bins;
    int out_of_range = 0;

    for (double value : data) {
        if (value < min_range || value >= max_range) {
            out_of_range++;
            continue;
        }

        int bin_idx = static_cast<int>((value - min_range) / bin_width);

        if (bin_idx < 0 || bin_idx >= num_bins) {
            out_of_range++;
            continue;
        }

        counts[bin_idx]++;
    }

    std::cout << out_of_range << " / " << data.size() << " out of range";

    std::vector<double> density(num_bins, 0.0);

    double normalization_factor = data.size() * bin_width;
    for (int i = 0; i < num_bins; ++i)
        density[i] = counts[i] / normalization_factor;

    return density;
}

std::vector<double> generate_positions(alfaStabilny& levyGenerator, int num_particles, double dt, int steps, double noise) {
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
            double drift = get_force(x) * dt;

            double jump = levyGenerator.losuj() * noise;

            x += drift + jump;
        }
        positions.push_back(x);
    }
    cout << "\rProgress: 100%" << std::endl;

    return positions;
}


vector<double> get_histogram() {
    std::random_device rd;
    std::mt19937 gen(rd());

    double alpha = 1.9;
    double beta = 0.0;
    double mi = 0.0;
    double sigma = 1.0;
    double length = 40.0;
    int num_intervals = 1000;

    alfaStabilny levyGenerator(alpha, beta, mi, sigma, &gen);

    int num_particles = 500000;
    double dt = 0.01;
    int steps = 100;
    double noise_scaling = std::pow(dt, 1.0 / alpha);


    auto positions = generate_positions(levyGenerator, num_particles, dt, steps, noise_scaling);
    auto histogram = generate_histogram(positions, num_intervals + 1, -length/2, length/2);

    return histogram;
}