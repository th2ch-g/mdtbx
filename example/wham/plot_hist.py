import matplotlib.pyplot as plt
import numpy as np


def plot_histograms(input_file="hist.xvg", output_file="hist_plot.png"):
    # xvgファイルを読み込む
    data = np.loadtxt(input_file, comments=["#", "@"])

    # 1列目が反応座標
    x = data[:, 0]

    # 2列目以降が各ウィンドウのヒストグラムデータ
    # data.shape[1] で列数を取得し、ループで全ウィンドウをプロット
    num_windows = data.shape[1] - 1

    fig, ax = plt.subplots(figsize=(10, 6))

    # カラーマップを使用（ウィンドウ数が多い場合に色を自動で変える）
    cmap = plt.get_cmap("jet")
    colors = [cmap(i) for i in np.linspace(0, 1, num_windows)]

    for i in range(num_windows):
        y = data[:, i + 1]
        # カウントが0だけのデータは描画しない場合などの判定を入れても良いが、
        # 通常はそのままプロットして問題ない
        ax.plot(x, y, color=colors[i], alpha=0.7, linewidth=1.5)

    # ラベル設定
    ax.set_xlabel("Reaction Coordinate (nm)", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("Umbrella Sampling Histograms", fontsize=16)

    # グリッドなどの体裁
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=12)

    # 保存
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

    # notitle版
    # ax.set_title("")
    # ax.set_xlabel("")
    # ax.set_ylabel("")
    # ax.tick_params(axis='both', labelbottom=False, labelleft=False)
    # plt.savefig(output_file.replace(".png", "_no_title.png"), dpi=300, transparent=True)

    plt.close(fig)


# 実行
if __name__ == "__main__":
    plot_histograms()
