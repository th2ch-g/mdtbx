import numpy as np
import pymbar
from pymbar import timeseries
import sys
from pathlib import Path

# =============================================================================
# 1. Metadata Generation
# =============================================================================
kconst = 100
subsample = True

path_list = sorted(
    Path(f"./").glob(f"rep*"),
    key=lambda x: float(x.name.split("_")[-2])
)

with open("metadata.dat", "w") as f:
    for trial in path_list:
        try:
            # target_distance = float(trial.name.split("_")[-2]) * 0.1
            target_dihedral = float(trial.name.split("_")[-2])
        except:
            continue

        # us1_pullx.xvg から us5... まで探す
        # ファイル名パターンがターゲットによって違う場合はここを調整してください
        # for us_trial in range(1, 5+1):
        # # for us_trial in [1]:
        #     filepath = f"{trial}/us{us_trial}_pullx.xvg"
        #     if Path(filepath).exists():
        #         f.write(f"{filepath} {target_dihedral} {kconst}\n")

        filepath = f"{trial}/us_dih1_pullx.xvg"
        if Path(filepath).exists():
            f.write(f"{filepath} {target_dihedral} {kconst}\n")

# =============================================================================
# Constants
# =============================================================================
kB = 8.314462618e-3
temperature = 310.0
beta = 1.0 / (kB * temperature)

# =============================================================================
# Read Data
# =============================================================================
files = []
r0_k = []
K_k = []
try:
    with open("metadata.dat", "r") as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 3:
                files.append(tokens[0])
                r0_k.append(float(tokens[1]))
                K_k.append(float(tokens[2]))
except FileNotFoundError:
    print("Error: metadata.dat could not be created.")
    sys.exit(1)

K = len(files)
if K == 0:
    print("Error: No data files found. Check your file paths and naming convention.")
    sys.exit(1)

r0_k = np.array(r0_k)
K_k = np.array(K_k)
N_max = 200000

N_k = np.zeros(K, dtype=int)
r_kn = np.zeros([K, N_max])
u_kn = np.zeros([K, N_max])

print(f"Reading {K} files...")

for k in range(K):
    filename = files[k]
    raw_data = []
    try:
        with open(filename, "r") as infile:
            for line in infile:
                if not line.startswith(("#", "@")):
                    parts = line.split()
                    if len(parts) >= 2:
                        val = float(parts[1])
                        # たまに極端な外れ値(0など)が入ることがあるのでフィルタリング
                        if val > 0.001:
                            raw_data.append(val)
        raw_data = np.array(raw_data)
    except FileNotFoundError:
        continue

    if len(raw_data) == 0:
        continue

    # --- 【修正1】サブサンプリングを無効化（全データ使用） ---
    # データ数が少ないため、間引き処理を行わずにすべてのサンプルを使います。
    if subsample:
        g = timeseries.statistical_inefficiency(raw_data)
        indices = timeseries.subsample_correlated_data(raw_data, g=g)

        # # t0: 平衡化開始地点, g: 統計的非効率性, Neff: 実効サンプル数
        # t0, g, Neff = timeseries.detect_equilibration(raw_data)
        #
        # # 平衡化前のデータを捨てて、かつ間引いたインデックスを取得
        # data_equil = raw_data[t0:]
        # indices_equil = timeseries.subsample_correlated_data(data_equil, g=g)
        #
        # # 元の配列(raw_data)に対するインデックスに変換
        # indices = indices_equil + t0
        #
        # # 確認用ログ（どれくらい間引かれたか確認できます）
        # print(f"Window {k}: Total={len(raw_data)}, EquilStart={t0}, g={g:.2f}, Kept={len(indices)}")
    else:
        indices = np.arange(len(raw_data))
    # -----------------------------------------------------

    N_k[k] = len(indices)

    if N_k[k] > r_kn.shape[1]:
        new_size = max(N_k[k], r_kn.shape[1] * 2)
        r_kn = np.pad(r_kn, ((0,0), (0, new_size - r_kn.shape[1])), 'constant')
        u_kn = np.pad(u_kn, ((0,0), (0, new_size - u_kn.shape[1])), 'constant')

    r_kn[k, 0:N_k[k]] = raw_data[indices]

    if k % 10 == 0:
        print(f"Processed window {k}/{K} (Samples: {N_k[k]})")

# 配列を実サイズに切り詰め
N_max = np.max(N_k)
if N_max == 0:
    print("Error: No valid data samples loaded.")
    sys.exit(1)

r_kn = r_kn[:, :N_max]
u_kn = u_kn[:, :N_max]

# =============================================================================
# Auto-detect Range & Prepare MBAR
# =============================================================================
valid_mask = r_kn > 0.001
valid_data = r_kn[valid_mask]
data_min = valid_data.min()
data_max = valid_data.max()

# --- 【修正2】範囲の余白（パディング）を極小にする ---
# 以前は 0.05 でしたが、空ビンを作らないよう 0.001 (1e-3) 程度にします
padding = 0.001
dist_min = data_min - padding
dist_max = data_max + padding

# ビン数の決定
range_width = dist_max - dist_min
# データ密度に応じてビン幅を調整（ここでは少し粗めの 0.03 nm 程度から試す）
target_bin_width = 0.03
nbins = int(range_width / target_bin_width)
nbins = max(5, min(nbins, 50)) # ビン数が少なすぎず多すぎないように制限

print(f"Auto-detected range: {dist_min:.4f} nm - {dist_max:.4f} nm")
print(f"Set nbins to: {nbins}")

print("Evaluating reduced potential energy matrix...")
u_kln = np.zeros([K, K, N_max])
for k in range(K):
    # k番目のシミュレーションのデータを取り出す
    r = r_kn[k, :N_k[k]]
    # 全ウィンドウ l でのエネルギーを計算
    # diff[l, n] = r[n] - r0_k[l]
    diff = r[np.newaxis, :] - r0_k[:, np.newaxis]
    u_kln[k, :, :N_k[k]] = u_kn[k, :N_k[k]] + beta * 0.5 * K_k[:, np.newaxis] * (diff**2)

# =============================================================================
# Run MBAR
# =============================================================================
print("Running MBAR (Robust mode)...")
# fes = pymbar.FES(u_kln, N_k, verbose=True, mbar_options={'solver_protocol': 'robust'})
fes = pymbar.FES(u_kln, N_k, verbose=True)

bin_edges = np.linspace(dist_min, dist_max, nbins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

r_n = pymbar.utils.kn_to_n(r_kn, N_k=N_k)
u_n = pymbar.utils.kn_to_n(u_kn, N_k=N_k)

histogram_parameters = {"bin_edges": bin_edges}
fes.generate_fes(u_n, r_n, fes_type="histogram", histogram_parameters=histogram_parameters)
# fes.generate_fes(u_kn, r_n, fes_type="histogram", histogram_parameters=histogram_parameters, n_bootstraps=100)

print("Computing FES...")
# 念のためエラー計算なしでトライ
try:
    results = fes.get_fes(bin_centers, reference_point="from-lowest", uncertainty_method="analytical")
    # results = fes.get_fes(bin_centers, reference_point="from-lowest", uncertainty_method="bootstrap")
    # results = fes.get_fes(bin_centers, reference_point="from-lowest", uncertainty_method="bootstrap", n_bootstraps=100, bootstrap_solver_protocol="robust")
except Exception as e:
    print(f"Analytical error calculation failed: {e}")
    print("Retrying without uncertainty calculation...")
    # results = fes.get_fes(bin_centers, reference_point="from-lowest", uncertainty_method=None)
    exit(1)

f_i = results["f_i"]
# 誤差が計算できなかった場合はゼロ埋め
df_i = results.get("df_i", np.zeros_like(f_i))

print("\nFree Energy Profile (Histogram) [unit: kT]")
print(f"{'Dist(nm)':>10s} {'PMF(kT)':>10s} {'Error':>10s}")
out_data = []
for i in range(nbins):
    # nan や inf を除外して出力
    if np.isfinite(f_i[i]):
        print(f"{bin_centers[i]:10.4f} {f_i[i]:10.4f} {df_i[i]:10.4f}")
        out_data.append([bin_centers[i], f_i[i], df_i[i]])

np.savetxt("pmf_output.dat", out_data, header="Dihedral(degree) PMF(kT) Error(kT)")
print("\nDone! Saved to pmf_output.dat")

