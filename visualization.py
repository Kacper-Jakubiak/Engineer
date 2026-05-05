import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            values = [float(x) for x in line.strip().split(';') if x]
            data.append(values)
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


def main():
    # Load data
    filename = "alfa1.txt"
    data = read_data(filename)

    plot_data_with_slider(data)

if __name__ == "__main__":
  main()

