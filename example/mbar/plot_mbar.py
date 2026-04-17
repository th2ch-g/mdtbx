#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt
import numpy as np

GAS_CONSTANT = 0.0019872041  # kcal/(mol·K)


def main():
    parser = argparse.ArgumentParser(description="Plot PMF profile from MBAR output")
    parser.add_argument("input_file", help="PMF data file (e.g. pmf_output.dat)")
    parser.add_argument(
        "--output", default="pmf_profile.png", help="Output image filename"
    )
    parser.add_argument(
        "--temperature", type=float, default=310.0, help="Temperature [K]"
    )
    parser.add_argument("--ymax", type=float, default=None, help="Y-axis upper limit")
    args = parser.parse_args()

    kT_to_kcal = GAS_CONSTANT * args.temperature
    print(f"Temperature: {args.temperature} K  |  1 kT = {kT_to_kcal:.4f} kcal/mol")

    try:
        data = np.loadtxt(args.input_file)
    except FileNotFoundError:
        print(f"Error: {args.input_file} not found.")
        raise SystemExit(1)

    r = data[:, 0]
    pmf_kcal = data[:, 1] * kT_to_kcal
    error_kcal = data[:, 2] * kT_to_kcal

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.fill_between(
        r,
        pmf_kcal - error_kcal,
        pmf_kcal + error_kcal,
        color="#1f77b4",
        alpha=0.3,
        linewidth=0,
        label="Standard Error",
    )
    ax.plot(r, pmf_kcal, color="#1f77b4", linewidth=2.5, label="PMF")
    ax.set_title(
        f"Potential of Mean Force (T={int(args.temperature)}K)",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_xlabel("Reaction Coordinate", fontsize=14)
    ax.set_ylabel("Free Energy (kcal/mol)", fontsize=14)
    if args.ymax is not None:
        ax.set_ylim(top=args.ymax)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=12, direction="in")
    ax.legend(fontsize=12, loc="best", frameon=True, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
