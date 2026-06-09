import matplotlib.pyplot as plt
import numpy as np

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


def compare_distributions(distribution1, distribution2):
    if distribution1.shape != distribution2.shape:
        raise ValueError(f"Shapes do not match: {distribution1.shape} != {distribution2.shape}")
    mae = np.mean(np.abs(distribution1 - distribution2))
    mse = np.mean((distribution1 - distribution2) ** 2)
    cdf_f = np.cumsum(distribution1)
    cdf_calc = np.cumsum(distribution2)
    ks_stat = np.max(np.abs(cdf_f - cdf_calc))
    print()
    print(f"MAE: {mae:.6f}")
    print(f"MSE: {mse:.6f}")
    print(f"KS Statistic: {ks_stat:.6f}")


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
    # print(f"{sum(f) = }")



    print(f"{dx = }")
    print(f"{1/dx = }")

    x_min, x_max = - L / 2, L / 2
    xs = np.linspace(x_min, x_max, points)

    compare_distributions(np.array(calculated), np.array(histogram))

    # plt.plot(xs, f, label="Inverse FFT")
    plt.plot(xs, histogram, label="Histogram")
    plt.plot(xs, calculated, label="Calculated")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
