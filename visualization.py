import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from scipy.stats import norm


def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = [float(x) for x in line.strip().split('\t') if x]
            data.append(values)
    return data


def normal_distribution(length, I, sigma):
    xs = np.linspace(0, length, I)
    ys = norm.pdf(xs, loc=length / 2, scale=sigma)
    return ys


def normal_distribution_over_time(length, I, dt, J, f0, K=1.0):
    data = [f0]
    normal = sum(f0) / 100
    for j in range(1, J):
        data.append(normal * normal_distribution(length, I, np.sqrt(K * dt * j)))
    return data


import numpy as np


def plot_time_snapshot(ax, data, time_index=0):
    y = data[time_index]
    x = list(range(len(y)))

    ax.clear()
    ax.plot(x, y)
    ax.set_title(f"Time step: {time_index}")
    ax.set_xlabel("Index")
    ax.set_ylabel("Value")


def plot_sum(ax, data):
    arr = np.array(data)

    start_sum = np.sum(arr[0])
    sums = np.sum(arr, axis=1) / start_sum

    ax.clear()
    ax.plot(sums)
    ax.set_ylabel("Sum")
    ax.set_title("Sum over Time")
    ax.set_ylim(-0.03, 1.03)


def plot_symmetry_mse(ax, data):
    arr = np.array(data)

    reversed_arr = arr[:, ::-1]
    mse_per_timestep = np.mean((arr - reversed_arr) ** 2, axis=1)

    ax.clear()
    ax.plot(mse_per_timestep)
    ax.set_xlabel("Time step")
    ax.set_ylabel("MSE")
    ax.set_title("Symmetry MSE")


def plot_data_with_slider(data):
    time_index = 0

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    # Initial plot using your reusable function
    plot_time_snapshot(ax, data, time_index)

    # Slider setup
    ax_slider = plt.axes((0.2, 0.1, 0.6, 0.03))
    slider = Slider(
        ax_slider,
        'Time',
        0,
        len(data) - 1,
        valinit=time_index,
        valstep=1
    )

    # Update function
    def update(val):
        t = int(slider.val)
        plot_time_snapshot(ax, data, t)
        fig.canvas.draw_idle()

    slider.on_changed(update)

    plt.show()


def dashboard(data):
    arr = np.array(data)
    time_index = 0

    # --- Precompute metrics ---
    start_sum = np.sum(arr[0])
    sums = np.sum(arr, axis=1) / start_sum

    reversed_arr = arr[:, ::-1]
    mse = np.mean((arr - reversed_arr) ** 2, axis=1)

    t_values = np.arange(len(data))

    # --- Layout (snapshot bigger) ---
    fig, axs = plt.subplot_mosaic([["snap", "sum"], ["snap", "mse"]])
    ax_snap, ax_sum, ax_mse = axs["snap"], axs["sum"], axs["mse"]
    plt.subplots_adjust(bottom=0.2)

    # --- Initial plots ---
    plot_time_snapshot(ax_snap, data, time_index)

    plot_sum(ax_sum, data)
    plot_symmetry_mse(ax_mse, data)

    # --- Add moving markers ---
    sum_dot, = ax_sum.plot(time_index, sums[time_index], 'o')
    mse_dot, = ax_mse.plot(time_index, mse[time_index], 'o')

    # --- Slider ---
    ax_slider = plt.axes((0.2, 0.08, 0.6, 0.03))
    slider = Slider(
        ax_slider,
        'Time',
        0,
        len(data) - 1,
        valinit=time_index,
        valstep=1
    )

    # --- Update ---
    def update(val):
        t = int(slider.val)

        # Update snapshot
        plot_time_snapshot(ax_snap, data, t)

        # Update dots
        sum_dot.set_data([t], [sums[t]])
        mse_dot.set_data([t], [mse[t]])

        fig.canvas.draw_idle()

    slider.on_changed(update)

    plt.show()


def main():
    # Load data
    filename = "history.txt"  # result_1.90.txt"
    # filename = input("Enter the filename: ")
    data = read_data(filename)
    # K = sum(data[0])
    # data2 = normal_distribution_over_time(40, 4001, 0.001, 1001, data[0])

    dashboard(data)
    # dashboard(data2)


if __name__ == "__main__":
    main()
