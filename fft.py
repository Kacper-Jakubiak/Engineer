import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, cauchy
from scipy.special import rel_entr
from visualization import read_data

def function(k, sigma, alpha, beta, mi):
    if alpha == 1:
        half = -2/np.pi * np.log(np.abs(k))
    else:
        half = np.tan(np.pi*alpha/2)

    first = -1j * k * mi

    middle = np.abs(sigma * k)**alpha

    return np.exp(first - middle * (1 + 1j * beta * np.sign(k) * half))

def alg(N, L, sigma, alpha, beta, mi):
    k = 2*np.pi * np.fft.fftfreq(N, d=L / N)

    g = function(k, sigma, alpha, beta, mi)

    f = np.fft.ifft(g)
    f = np.real(f)

    f = np.fft.fftshift(f)
    x_plus = np.linspace(-L/2, L/2, N+1)#, endpoint=False)
    f_plus = np.append(f, f[0])

    return x_plus, f_plus

def compare_distributions(calculated, theoretical):
    if theoretical.shape != calculated.shape:
        raise ValueError(f"Shapes do not match: {theoretical.shape} != {calculated.shape}")
    # theoretical = theoretical / np.sum(theoretical)
    # calculated = calculated / np.sum(calculated)
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
    filename = "K_1.00.txt"#result_1.90.txt"
    L = 40
    K = 1.0
    dt = 0.001
    J = 1000
    sigma = K * dt * J
    beta = 0.0
    alpha = 1.9
    mi = 0.0

    data = read_data(filename)
    N = len(data[-1]) - 1
    print(f"{N = }")


    calculated = np.array(data[-1])
    calculated /= (N + 1)

    x, f = alg(N, L, sigma, alpha, beta, mi)

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