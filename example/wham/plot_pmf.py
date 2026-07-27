#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def plot_pmf(input_file="profile.xvg", output_file="pmf_plot.png"):
    # Read the xvg file (lines starting with # or @ are treated as comments)
    data = np.loadtxt(input_file, comments=["#", "@"])

    # Column 1 is the reaction coordinate (nm), column 2 is the free energy
    x = data[:, 0]
    y = data[:, 1]
    y = y - np.min(y)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot
    ax.plot(x, y, color="black", linewidth=2, label="PMF")

    # Labels
    ax.set_xlabel("Reaction Coordinate (nm)", fontsize=14)
    ax.set_ylabel("PMF (kcal/mol)", fontsize=14)
    ax.set_title("Potential of Mean Force", fontsize=16)

    # Grid and other cosmetics
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=12)

    # Save
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


# Run
if __name__ == "__main__":
    plot_pmf()
