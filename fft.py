import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, cauchy
from scipy.special import rel_entr
from visualization import read_data

def alg(N, L, sigma, alpha):
    dx = L / N

    x = np.arange(N) * dx
    k = 2*np.pi * np.fft.fftfreq(N, d=dx)

    g = np.exp(-np.abs(sigma * k)**alpha)

    dk = 2*np.pi / L
    f = np.fft.ifft(g) * N * dk / (2*np.pi)
    f = np.real(f)

    f = np.fft.fftshift(f)
    x = x - L/2

    return x, f

def compare_distributions(calculated, theoretical):
    mae = np.mean(np.abs(theoretical - calculated))
    mse = np.mean((theoretical - calculated)**2)
    kl_div = np.sum(rel_entr(calculated, theoretical))
    cdf_f = np.cumsum(theoretical)
    cdf_calc = np.cumsum(calculated)
    ks_stat = np.max(np.abs(cdf_f - cdf_calc))
    print(f"MAE: {mae:.6f}")
    print(f"MSE: {mse:.6f}")
    print(f"KL Divergence: {kl_div:.6f}")
    print(f"KS Statistic: {ks_stat:.6f}")

def main():
    filename = "output.txt"
    L = 100
    sigma = 1.0 * 0.001 * 2500
    alpha = 0.9
    data = read_data(filename)
    N = len(data[0])
    print(N)


    calculated = np.array(data[-1])
    calculated /= N

    x, f = alg(N, L, sigma, alpha)
    f *= L / N

    print(np.sum(f))
    print(np.sum(calculated))
    compare_distributions(calculated, f)

    # gaussian = norm.pdf(x, loc=0, scale=np.sqrt(2))
    # cauchy_dist = cauchy.pdf(x, loc=0, scale=1)

    # plt.plot(x, gaussian, label="Gaussian (SciPy)")
    # plt.plot(x, cauchy_dist, label="Cauchy (SciPy)")
    plt.plot(x, f, label="Inverse FFT")  
    plt.plot(x, calculated, label="Calculated")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()