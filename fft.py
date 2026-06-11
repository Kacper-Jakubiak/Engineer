import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import jensenshannon


def read_line(filename) -> list[float]:
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = [float(x) for x in line.strip().split('\t') if x]
            break
    return values


def get_length(filename) -> float:
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#length'):
                return float(line.strip().split(' ')[-1])
    raise ValueError("Length not found in file")


def function(k, sigma, alpha, beta, mi):
    if alpha == 1:
        theta = -2 / np.pi * np.log(np.abs(k))
    else:
        theta = np.tan(np.pi * alpha / 2)

    first = -1j * k * mi

    middle = np.abs(sigma ** (1.0 / alpha) * k) ** alpha

    last = 1 - -1j * beta * np.sign(k) * theta

    return np.exp(first - middle * last)


def alg(points, L, sigma, alpha, beta, mi):
    N = points - 1
    k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)

    g = function(k, sigma, alpha, beta, mi)

    f = np.fft.ifft(g)
    f = np.real(f)

    f = np.fft.fftshift(f)
    f_plus = np.append(f, f[0])

    return f_plus.tolist()


def compare_distributions(dist1, dist2, dx):
    dist1 = np.asarray(dist1)
    dist2 = np.asarray(dist2)
    if dist1.shape != dist2.shape:
        raise ValueError(f"Shapes do not match: {dist1.shape} != {dist2.shape}")

    mae = np.mean(np.abs(dist1 - dist2))
    mse = np.mean((dist1 - dist2) ** 2)

    cdf1 = np.cumsum(dist1)
    cdf2 = np.cumsum(dist2)
    ks_stat = np.max(np.abs(cdf1 - cdf2))

    tvd = 0.5 * np.sum(np.abs(dist1 - dist2))

    ovl = np.sum(np.minimum(dist1, dist2)) * dx

    hellinger = np.sqrt(0.5 * np.sum((np.sqrt(dist1) - np.sqrt(dist2)) ** 2))

    print(f"MAE: {mae:.6}, MSE: {mse:.6}")
    print(f"KS Statistic: {ks_stat:.6}")
    print(f"Total Variation Distance: {tvd:.6}")
    print(f"Hellinger Distance: {hellinger:.6}")
    print(f"Overlap Coefficient: {ovl:.6}")


def main():
    filename = "result.csv"
    histname = "histogram.csv"
    L = get_length(filename)
    print(f"{L = }")

    sigma = 1.0
    alpha = 1.9
    beta = 0.0
    mi = 0.0

    time = 1.0

    calculated = read_line(filename)
    histogram = read_line(histname)

    points = len(calculated)
    dx = L / points

    # f = alg(points, L, sigma * time, alpha, beta, mi * time)
    # f = [f_i * points for f_i in f]

    # print(f"{len(calculated) = }")
    # print(f"{len(histogram) = }")
    print(f"{sum(calculated) = }")
    print(f"{sum(histogram) = }")
    print(f"{sum(calculated) * dx = }")
    print(f"{sum(histogram) * dx = }")
    # print(f"{sum(f) = }")

    print(f"{dx = }")
    print(f"{1/dx = }")

    x_min, x_max = - L / 2, L / 2
    xs = np.linspace(x_min, x_max, points)

    compare_distributions(calculated, histogram, dx)

    # plt.plot(xs, f, label="Inverse FFT")
    plt.plot(xs, histogram, label="Histogram")
    plt.plot(xs, calculated, label="Calculated")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
