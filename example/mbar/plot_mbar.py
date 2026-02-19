import numpy as np
import matplotlib.pyplot as plt
import sys

# ==========================================
# 設定
# ==========================================
input_file = sys.argv[1]
output_image = "pmf_profile.png"

# 温度設定 (計算時と同じ温度を指定してください)
temperature = 310.0  # Kelvin

# 単位変換係数の計算 (kT -> kcal/mol)
# ガス定数 R = 0.0019872 kcal/(mol·K)
gas_constant = 0.0019872041
kT_to_kcal = gas_constant * temperature

print(f"Temperature: {temperature} K")
print(f"Conversion factor (1 kT): {kT_to_kcal:.4f} kcal/mol")

# ==========================================
# データの読み込み
# ==========================================
try:
    data = np.loadtxt(input_file)
    r = data[:, 0]  # 1列目: 距離 (nm)
    pmf_kT = data[:, 1]  # 2列目: PMF (kT)
    error_kT = data[:, 2]  # 3列目: 誤差 (kT)
except FileNotFoundError:
    print(f"Error: {input_file} が見つかりません。")
    exit()

# 単位変換を実行
pmf_kcal = pmf_kT * kT_to_kcal
error_kcal = error_kT * kT_to_kcal

# ==========================================
# プロットの作成
# ==========================================
fig, ax = plt.subplots(figsize=(8, 6))

# 1. 誤差範囲を塗りつぶし (Shaded Error Bar)
# kcal/mol なので少し線や色を濃いめに見やすく設定
ax.fill_between(
    r,
    pmf_kcal - error_kcal,
    pmf_kcal + error_kcal,
    color="#1f77b4",
    alpha=0.3,
    linewidth=0,
    label="Standard Error",
)

# 2. PMFのメインライン
ax.plot(r, pmf_kcal, color="#1f77b4", linewidth=2.5, label="PMF")

# ==========================================
# 装飾
# ==========================================
ax.set_title(
    f"Potential of Mean Force (T={int(temperature)}K)", fontsize=16, fontweight="bold"
)
ax.set_xlabel("Dihedral (nm)", fontsize=14)
ax.set_ylabel("Free Energy (kcal/mol)", fontsize=14)

plt.ylim(0, 9)

# # ゼロライン（基準線）を引く
# ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

# グリッド線
ax.grid(True, linestyle=":", alpha=0.6)

# 軸の文字サイズと目盛りの向き
ax.tick_params(axis="both", which="major", labelsize=12, direction="in")

# 凡例
ax.legend(fontsize=12, loc="best", frameon=True, framealpha=0.9)

# 余白の調整
plt.tight_layout()

# ==========================================
# 保存
# ==========================================
plt.savefig(output_image, dpi=300)
print(f"グラフを保存しました: {output_image}")

ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="both", labelbottom=False, labelleft=False)
legend = ax.get_legend()
if legend:
    legend.remove()
plt.savefig(output_image.replace(".png", "_no_title.png"), dpi=300, transparent=True)
