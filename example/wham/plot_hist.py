#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def plot_histograms(input_file="hist.xvg", output_file="hist_plot.png"):
    # Read the xvg file
    data = np.loadtxt(input_file, comments=["#", "@"])

    # Column 1 is the reaction coordinate
    x = data[:, 0]

    # Columns 2 onwards hold the histogram of each window
    # data.shape[1] gives the column count; loop over it to plot every window
    num_windows = data.shape[1] - 1

    fig, ax = plt.subplots(figsize=(10, 6))

    # Use a colour map (colours vary automatically when there are many windows)
    cmap = plt.get_cmap("jet")
    colors = [cmap(i) for i in np.linspace(0, 1, num_windows)]

    for i in range(num_windows):
        y = data[:, i + 1]
        # A check to skip data whose counts are all zero could be added here,
        # but plotting it as is normally causes no problem
        ax.plot(x, y, color=colors[i], alpha=0.7, linewidth=1.5)

    # Labels
    ax.set_xlabel("Reaction Coordinate (nm)", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("Umbrella Sampling Histograms", fontsize=16)

    # Grid and other cosmetics
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=12)

    # Save
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


# Run
if __name__ == "__main__":
    plot_histograms()
