import matplotlib.pyplot as plt
import numpy as np

def plot_pmf(input_file='profile.xvg', output_file='pmf_plot.png'):
    # xvgファイルを読み込む（#と@で始まる行はコメントとして無視）
    data = np.loadtxt(input_file, comments=['#', '@'])

    # 1列目が反応座標(nm)、2列目が自由エネルギー
    x = data[:, 0]
    y = data[:, 1]
    y = y - np.min(y)

    fig, ax = plt.subplots(figsize=(8, 6))

    # プロット
    ax.plot(x, y, color='black', linewidth=2, label='PMF')

    # ラベル設定
    ax.set_xlabel('Reaction Coordinate (nm)', fontsize=14)
    ax.set_ylabel('PMF (kcal/mol)', fontsize=14)
    ax.set_title('Potential of Mean Force', fontsize=16)

    # グリッドなどの体裁
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # 保存
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

    # ax.set_title("")
    # ax.set_xlabel("")
    # ax.set_ylabel("")
    # ax.tick_params(axis='both', labelbottom=False, labelleft=False)
    # plt.savefig(output_file.replace(".png", "_no_title.png"), dpi=300, transparent=True)


    plt.close(fig) # メモリ解放のためclose

# 実行
if __name__ == "__main__":
    plot_pmf()
