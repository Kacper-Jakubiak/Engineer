import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, cauchy

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

def main():
    N = 2**12
    L = 40
    sigma = 1.0
    alpha = 1.0

    x, f = alg(N, L, sigma, alpha)

    gaussian = norm.pdf(x, loc=0, scale=np.sqrt(2))
    cauchy_dist = cauchy.pdf(x, loc=0, scale=1)

    plt.plot(x, gaussian, label="Gaussian (SciPy)")
    plt.plot(x, cauchy_dist, label="Cauchy (SciPy)")
    plt.plot(x, f, label="Inverse FFT")  
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()