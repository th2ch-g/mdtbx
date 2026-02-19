import numpy as np
import sys
import matplotlib.pyplot as plt

# =============================================================================
# Constants
# =============================================================================
# 前回のスクリプトと同じ温度設定
temperature = 310.0

# 単位変換係数
# Gas constant R in kcal/(mol*K)
R_kcal = 1.987204259e-3
RT_kcal = R_kcal * temperature

# =============================================================================
# Arguments & Usage
# =============================================================================
if len(sys.argv) < 3:
    print("Usage: python calc_delta_g.py <dist_A_nm> <dist_B_nm> [input_file]")
    print("Example: python calc_delta_g.py 1.2 1.7")
    sys.exit(1)

dist_A_target = float(sys.argv[1])
dist_B_target = float(sys.argv[2])
input_file = sys.argv[3] if len(sys.argv) > 3 else "pmf_output.dat"

# =============================================================================
# Load Data
# =============================================================================
print(f"Loading {input_file} ...")
try:
    # 読み込み (headerの#は自動的に無視される)
    data = np.loadtxt(input_file)
except OSError:
    print(f"Error: File {input_file} not found.")
    sys.exit(1)

dist = data[:, 0]  # Distance (nm)
pmf_kt = data[:, 1]  # PMF (kT)
err_kt = data[:, 2]  # Error (kT)

# =============================================================================
# Find Nearest Points
# =============================================================================
# 指定された距離に最も近いビンのインデックスを探す
idx_A = (np.abs(dist - dist_A_target)).argmin()
idx_B = (np.abs(dist - dist_B_target)).argmin()

dist_A_actual = dist[idx_A]
pmf_A_kt = pmf_kt[idx_A]
err_A_kt = err_kt[idx_A]

dist_B_actual = dist[idx_B]
pmf_B_kt = pmf_kt[idx_B]
err_B_kt = err_kt[idx_B]

# =============================================================================
# Calculate Delta G
# =============================================================================
# Delta G (kT) = PMF(B) - PMF(A)
# つまり、AからBへ移行する際のエネルギー変化
dG_kt = pmf_B_kt - pmf_A_kt

# 誤差伝播 (二乗和の平方根)
error_dG_kt = np.sqrt(err_A_kt**2 + err_B_kt**2)

# Convert to kcal/mol
pmf_kcal = pmf_kt * RT_kcal
dG_kcal = dG_kt * RT_kcal
error_dG_kcal = error_dG_kt * RT_kcal

# =============================================================================
# Output Results
# =============================================================================
print("-" * 60)
print(f"Temperature: {temperature} K")
print(f"RT factor  : {RT_kcal:.4f} kcal/mol")
print("-" * 60)
print(f"{'Point':<10} | {'Target(nm)':<10} | {'Actual(nm)':<10} | {'PMF(kT)':<10} | {'PMF(kcal/mol)':<15}")
print("-" * 60)
print(f"{'Start (A)':<10} | {dist_A_target:<10.4f} | {dist_A_actual:<10.4f} | {pmf_A_kt:<10.4f} | {pmf_A_kt * RT_kcal:<15.4f}")
print(f"{'End   (B)':<10} | {dist_B_target:<10.4f} | {dist_B_actual:<10.4f} | {pmf_B_kt:<10.4f} | {pmf_B_kt * RT_kcal:<15.4f}")
print("-" * 60)
print(f"Delta G (A -> B): {dG_kcal:.4f} +/- {error_dG_kcal:.4f} kcal/mol")
print("-" * 60)

# =============================================================================
# Plotting
# =============================================================================
plt.figure(figsize=(8, 5))
plt.errorbar(dist, pmf_kcal, yerr=err_kt*RT_kcal, fmt='-', color='black', ecolor='lightgray', label='PMF Profile')

# ポイントのハイライト
plt.scatter([dist_A_actual], [pmf_A_kt * RT_kcal], color='blue', s=100, zorder=5, label='Start (A)')
plt.scatter([dist_B_actual], [pmf_B_kt * RT_kcal], color='red', s=100, zorder=5, label='End (B)')

plt.title(f"PMF Profile with Selected Points\nDelta G = {dG_kcal:.2f} kcal/mol")
plt.xlabel("Distance (nm)")
plt.ylabel("PMF (kcal/mol)")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# 画像保存
plot_filename = "delta_g_check.png"
plt.savefig(plot_filename, dpi=150)
print(f"Plot saved to: {plot_filename}")
