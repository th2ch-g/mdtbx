#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def plot_bootstrap_pmf(result_file="bsResult.xvg", output_file="pmf_bootstrap.png"):
    # --- Read bsResult.xvg (mean and standard deviation) ---
    # column 1: reaction coordinate, column 2: mean PMF, column 3: standard deviation
    data = np.loadtxt(result_file, comments=["#", "@"])

    x = data[:, 0]
    y_mean = data[:, 1]
    y_std = data[:, 2]

    # --- Shift the energy minimum to zero ---
    # Find the minimum of the mean profile and shift the whole curve by it
    min_idx = np.argmin(y_mean)
    min_val = y_mean[min_idx]

    y_mean_shifted = y_mean - min_val

    # --- Build the plot ---
    fig, ax = plt.subplots(figsize=(8, 6))

    # 1. Draw the standard deviation as a shaded band (Mean +/- Std)
    # upper and lower bounds
    upper_bound = y_mean_shifted + y_std
    lower_bound = y_mean_shifted - y_std

    ax.fill_between(
        x, lower_bound, upper_bound, color="gray", alpha=0.4, label="Standard Deviation"
    )

    # 2. Draw the mean PMF as a line
    ax.plot(x, y_mean_shifted, color="black", linewidth=2, label="Mean PMF")

    # --- Tidy up the appearance ---
    ax.set_xlabel("Reaction Coordinate (nm)", fontsize=14)
    ax.set_ylabel("PMF (kcal/mol)", fontsize=14)
    ax.set_title("PMF with Bootstrap Error (n=200)", fontsize=16)

    ax.grid(True, linestyle=":", alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.legend(loc="best", fontsize=12)

    # Save
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)
    print(f"Saved plot to {output_file}")


if __name__ == "__main__":
    plot_bootstrap_pmf()
