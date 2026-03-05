import sys
from pathlib import Path

import numpy as np
import pymbar
from pymbar import timeseries

# =============================================================================
# Configuration
# =============================================================================
KCONST = 100  # Umbrella spring constant [kJ/mol/unit^2]
SUBSAMPLE = True  # Statistical subsampling to reduce autocorrelation
TEMPERATURE = 310.0  # [K]
TARGET_BIN_WIDTH = 0.03
NBINS_MIN = 5
NBINS_MAX = 50
PADDING = 0.001  # Range padding to avoid empty edge bins
N_MAX_INIT = 200000  # Initial array size for trajectory data


def read_metadata(metadata_file: str) -> tuple[list, list, list]:
    files, r0_k, K_k = [], [], []
    with open(metadata_file) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 3:
                files.append(tokens[0])
                r0_k.append(float(tokens[1]))
                K_k.append(float(tokens[2]))
    return files, r0_k, K_k


def build_metadata(path_list: list[Path], kconst: float) -> str:
    lines = []
    for trial in path_list:
        try:
            target_dihedral = float(trial.name.split("_")[-2])
        except (ValueError, IndexError):
            continue
        filepath = f"{trial}/us_dih1_pullx.xvg"
        if Path(filepath).exists():
            lines.append(f"{filepath} {target_dihedral} {kconst}")
    content = "\n".join(lines)
    with open("metadata.dat", "w") as f:
        f.write(content)
    return content


def load_trajectories(files: list, r0_k: list, K_k: list, beta: float, subsample: bool):
    K = len(files)
    r0_k = np.array(r0_k)
    K_k = np.array(K_k)

    N_k = np.zeros(K, dtype=int)
    r_kn = np.zeros([K, N_MAX_INIT])
    u_kn = np.zeros([K, N_MAX_INIT])

    print(f"Reading {K} files...")
    for k, filename in enumerate(files):
        raw_data = []
        try:
            with open(filename) as infile:
                for line in infile:
                    if not line.startswith(("#", "@")):
                        parts = line.split()
                        if len(parts) >= 2:
                            val = float(parts[1])
                            if val > 0.001:
                                raw_data.append(val)
        except FileNotFoundError:
            continue

        raw_data = np.array(raw_data)
        if len(raw_data) == 0:
            continue

        if subsample:
            g = timeseries.statistical_inefficiency(raw_data)
            indices = timeseries.subsample_correlated_data(raw_data, g=g)
        else:
            indices = np.arange(len(raw_data))

        N_k[k] = len(indices)

        if N_k[k] > r_kn.shape[1]:
            new_size = max(N_k[k], r_kn.shape[1] * 2)
            r_kn = np.pad(r_kn, ((0, 0), (0, new_size - r_kn.shape[1])), "constant")
            u_kn = np.pad(u_kn, ((0, 0), (0, new_size - u_kn.shape[1])), "constant")

        r_kn[k, : N_k[k]] = raw_data[indices]

        if k % 10 == 0:
            print(f"  Window {k}/{K} (samples: {N_k[k]})")

    N_max = np.max(N_k)
    if N_max == 0:
        raise RuntimeError("No valid data samples loaded.")

    return r_kn[:, :N_max], u_kn[:, :N_max], N_k, r0_k, K_k


def run_mbar(r_kn, u_kn, N_k, r0_k, K_k, beta):
    valid_data = r_kn[r_kn > 0.001]
    dist_min = valid_data.min() - PADDING
    dist_max = valid_data.max() + PADDING

    nbins = int((dist_max - dist_min) / TARGET_BIN_WIDTH)
    nbins = max(NBINS_MIN, min(nbins, NBINS_MAX))

    print(f"Range: {dist_min:.4f} - {dist_max:.4f}, bins: {nbins}")

    N_max = r_kn.shape[1]
    K = len(N_k)
    u_kln = np.zeros([K, K, N_max])
    for k in range(K):
        r = r_kn[k, : N_k[k]]
        diff = r[np.newaxis, :] - r0_k[:, np.newaxis]
        u_kln[k, :, : N_k[k]] = u_kn[k, : N_k[k]] + beta * 0.5 * K_k[:, np.newaxis] * (
            diff**2
        )

    print("Running MBAR...")
    fes = pymbar.FES(u_kln, N_k, verbose=True)

    bin_edges = np.linspace(dist_min, dist_max, nbins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    r_n = pymbar.utils.kn_to_n(r_kn, N_k=N_k)
    u_n = pymbar.utils.kn_to_n(u_kn, N_k=N_k)
    fes.generate_fes(
        u_n, r_n, fes_type="histogram", histogram_parameters={"bin_edges": bin_edges}
    )

    results = fes.get_fes(
        bin_centers, reference_point="from-lowest", uncertainty_method="analytical"
    )
    f_i = results["f_i"]
    df_i = results.get("df_i", np.zeros_like(f_i))

    return bin_centers, f_i, df_i


def main():
    kB = 8.314462618e-3  # kJ/(mol·K)
    beta = 1.0 / (kB * TEMPERATURE)

    # Build metadata from directory listing
    path_list = sorted(
        Path("./").glob("rep*"), key=lambda x: float(x.name.split("_")[-2])
    )
    if path_list:
        build_metadata(path_list, KCONST)

    try:
        files, r0_k, K_k = read_metadata("metadata.dat")
    except FileNotFoundError:
        print("Error: metadata.dat not found.")
        sys.exit(1)

    if not files:
        print("Error: No data files found. Check your naming convention.")
        sys.exit(1)

    r_kn, u_kn, N_k, r0_k, K_k = load_trajectories(files, r0_k, K_k, beta, SUBSAMPLE)
    bin_centers, f_i, df_i = run_mbar(r_kn, u_kn, N_k, r0_k, K_k, beta)

    print(f"\n{'Coord':>10s} {'PMF(kT)':>10s} {'Error':>10s}")
    out_data = []
    for i in range(len(bin_centers)):
        if np.isfinite(f_i[i]):
            print(f"{bin_centers[i]:10.4f} {f_i[i]:10.4f} {df_i[i]:10.4f}")
            out_data.append([bin_centers[i], f_i[i], df_i[i]])

    np.savetxt("pmf_output.dat", out_data, header="Coord PMF(kT) Error(kT)")
    print("\nDone! Saved to pmf_output.dat")


if __name__ == "__main__":
    main()
