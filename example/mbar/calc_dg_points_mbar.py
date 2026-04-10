import argparse

import matplotlib.pyplot as plt
import numpy as np

R_KCAL = 1.987204259e-3  # Gas constant [kcal/(mol*K)]


def calc_delta_g(dist, pmf_kt, err_kt, dist_A_target, dist_B_target, temperature):
    RT_kcal = R_KCAL * temperature

    idx_A = (np.abs(dist - dist_A_target)).argmin()
    idx_B = (np.abs(dist - dist_B_target)).argmin()

    dist_A_actual = dist[idx_A]
    pmf_A_kt = pmf_kt[idx_A]
    err_A_kt = err_kt[idx_A]

    dist_B_actual = dist[idx_B]
    pmf_B_kt = pmf_kt[idx_B]
    err_B_kt = err_kt[idx_B]

    dG_kt = pmf_B_kt - pmf_A_kt
    error_dG_kt = np.sqrt(err_A_kt**2 + err_B_kt**2)
    dG_kcal = dG_kt * RT_kcal
    error_dG_kcal = error_dG_kt * RT_kcal

    return (
        dist_A_actual,
        pmf_A_kt,
        dist_B_actual,
        pmf_B_kt,
        dG_kcal,
        error_dG_kcal,
        RT_kcal,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Calculate delta G between two points on a PMF profile"
    )
    parser.add_argument("dist_A", type=float, help="Start point distance (nm)")
    parser.add_argument("dist_B", type=float, help="End point distance (nm)")
    parser.add_argument(
        "input_file", nargs="?", default="pmf_output.dat", help="PMF data file"
    )
    parser.add_argument(
        "--temperature", type=float, default=310.0, help="Temperature [K]"
    )
    parser.add_argument(
        "--output", default="delta_g_check.png", help="Output plot filename"
    )
    args = parser.parse_args()

    try:
        data = np.loadtxt(args.input_file)
    except OSError:
        print(f"Error: File {args.input_file} not found.")
        raise SystemExit(1)

    dist = data[:, 0]
    pmf_kt = data[:, 1]
    err_kt = data[:, 2]

    dist_A, pmf_A_kt, dist_B, pmf_B_kt, dG_kcal, error_dG_kcal, RT_kcal = calc_delta_g(
        dist, pmf_kt, err_kt, args.dist_A, args.dist_B, args.temperature
    )
    pmf_kcal = pmf_kt * RT_kcal

    print(f"{'─' * 60}")
    print(f"Temperature: {args.temperature} K  |  RT: {RT_kcal:.4f} kcal/mol")
    print(f"{'─' * 60}")
    print(
        f"{'Point':<10} {'Target(nm)':<12} {'Actual(nm)':<12} "
        f"{'PMF(kT)':<10} {'PMF(kcal/mol)':<15}"
    )
    print(f"{'─' * 60}")
    print(
        f"{'Start (A)':<10} {args.dist_A:<12.4f} {dist_A:<12.4f} "
        f"{pmf_A_kt:<10.4f} {pmf_A_kt * RT_kcal:<15.4f}"
    )
    print(
        f"{'End   (B)':<10} {args.dist_B:<12.4f} {dist_B:<12.4f} "
        f"{pmf_B_kt:<10.4f} {pmf_B_kt * RT_kcal:<15.4f}"
    )
    print(f"{'─' * 60}")
    print(f"Delta G (A -> B): {dG_kcal:.4f} +/- {error_dG_kcal:.4f} kcal/mol")
    print(f"{'─' * 60}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(
        dist,
        pmf_kcal,
        yerr=err_kt * RT_kcal,
        fmt="-",
        color="black",
        ecolor="lightgray",
        label="PMF Profile",
    )
    ax.scatter(
        [dist_A], [pmf_A_kt * RT_kcal], color="blue", s=100, zorder=5, label="Start (A)"
    )
    ax.scatter(
        [dist_B], [pmf_B_kt * RT_kcal], color="red", s=100, zorder=5, label="End (B)"
    )
    ax.set_title(f"PMF Profile\nΔG = {dG_kcal:.2f} kcal/mol")
    ax.set_xlabel("Distance (nm)")
    ax.set_ylabel("PMF (kcal/mol)")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Plot saved to: {args.output}")


if __name__ == "__main__":
    main()
