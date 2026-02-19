import matplotlib.pyplot as plt
import numpy as np

def plot_bootstrap_pmf(result_file='bsResult.xvg', output_file='pmf_bootstrap.png'):
    # --- bsResult.xvg (平均と標準偏差) の読み込み ---
    # 1列目: 反応座標, 2列目: 平均PMF, 3列目: 標準偏差
    data = np.loadtxt(result_file, comments=['#', '@'])

    x = data[:, 0]
    y_mean = data[:, 1]
    y_std  = data[:, 2]

    # --- エネルギー最小値を0に合わせる補正 ---
    # 平均プロファイルの最小値を見つけ、全体をシフトする
    min_idx = np.argmin(y_mean)
    min_val = y_mean[min_idx]

    y_mean_shifted = y_mean - min_val

    # --- プロット作成 ---
    fig, ax = plt.subplots(figsize=(8, 6))

    # 1. 標準偏差を「帯（影）」として描画 (Mean ± Std)
    # 上限と下限
    upper_bound = y_mean_shifted + y_std
    lower_bound = y_mean_shifted - y_std

    ax.fill_between(x, lower_bound, upper_bound, color='gray', alpha=0.4, label='Standard Deviation')

    # 2. 平均PMFを線で描画
    ax.plot(x, y_mean_shifted, color='black', linewidth=2, label='Mean PMF')

    # --- 体裁を整える ---
    ax.set_xlabel('Reaction Coordinate (nm)', fontsize=14)
    ax.set_ylabel('PMF (kcal/mol)', fontsize=14)
    ax.set_title('PMF with Bootstrap Error (n=200)', fontsize=16)

    ax.grid(True, linestyle=':', alpha=0.6)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc='best', fontsize=12)

    # 保存
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

    # ax.set_title("")
    # ax.set_xlabel("")
    # ax.set_ylabel("")
    # ax.tick_params(axis='both', labelbottom=False, labelleft=False)
    # plt.savefig(output_file.replace(".png", "_no_title.png"), dpi=300, transparent=True)

    plt.close(fig)
    print(f"Saved plot to {output_file}")

if __name__ == "__main__":
    plot_bootstrap_pmf()
