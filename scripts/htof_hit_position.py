import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# ------------------------------
# グローバル設定（スタイルなど）
# ------------------------------
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
    'font.size': 28,
    'axes.linewidth': 1.0,
    'axes.grid': False,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'xtick.minor.size': 5,
    'ytick.minor.size': 5
})


def plot(arg_dict):
    # ------------------------------
    # ROOTファイル読み込み
    # ------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_file_path = os.path.join(script_dir, arg_dict["data"])
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    # ------------------------------
    # x, y データ結合（セグメントごとに結合）
    # ------------------------------
    x_all = np.array([])
    y_all = np.array([])
    for i in range(34):
        # x = -tree["x"][i]  # 左右反転
        # y = tree["y"][i]
        x = -tree["x_high"][i]  # 左右反転
        y = tree["y_high"][i]
        x_all = np.hstack([x_all, x])
        y_all = np.hstack([y_all, y])

    # ------------------------------
    # ヒストグラム用のbin定義（細かく）
    # ------------------------------
    n_split = 5  # 各セグメント内をこの数で細分割
    seg_edges = [-140, -70, 0, 70] + list(np.arange(140, 70 + (33 - 4) * 70 + 1, 70))

    # x方向：セグメントごとの細かいbin生成
    x_edges_fine = []
    for i in range(len(seg_edges) - 1):
        x_sub = np.linspace(seg_edges[i], seg_edges[i + 1], n_split + 1)[:-1]
        x_edges_fine.extend(x_sub)
    x_edges_fine.append(seg_edges[-1])

    # y方向：均等分割
    y_edges_fine = np.linspace(-400, 400, 201)

    # ------------------------------
    # 2Dヒストグラムプロット
    # ------------------------------
    fig, ax = plt.subplots(figsize=(8, 8))
    counts, xedges, yedges, im = ax.hist2d(
        x_all, y_all,
        bins=(x_edges_fine, y_edges_fine),
        cmap='viridis',
        norm=Normalize(),
        cmin=1,
        zorder=100
    )

    # ------------------------------
    # セグメント番号ラベル設定（2個ごと）
    # ------------------------------
    centers = [(seg_edges[i] + seg_edges[i+1]) / 2 for i in range(len(seg_edges) - 1)]
    labels = ["1", "", "4, 5"] + [str(i+1) if i % 2 == 0 else "" for i in range(5, 34)]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=90, fontsize = 22)

    # ------------------------------
    # セグメント境界線（縦線）
    # ------------------------------
    for xg in seg_edges:
        ax.axvline(x=xg, color="k", linestyle=":", linewidth=1.0)

    # ------------------------------
    # 一部範囲の塗りつぶし（背景色）
    # ------------------------------
    ax.fill_betweenx(
        y=[-500, 500],
        x1=seg_edges[7],
        x2=seg_edges[27],
        color="red",
        alpha=0.1,
        zorder=0
    )

    # ------------------------------
    # 軸ラベル・範囲設定
    # ------------------------------
    ax.set_title(arg_dict["title"])
    ax.set_xlabel("Segment")
    ax.set_ylabel("y [mm]")
    ax.set_ylim(-450, 450)

    # ------------------------------
    # 仕上げと表示
    # ------------------------------
    plt.subplots_adjust(left = 0.18, right=0.98, top=0.9, bottom = 0.15)
    img_save_path = os.path.join(
        script_dir, 
        "../results/img/htof/htof_hit_position_{}.png".format(os.path.splitext(os.path.basename(arg_dict["data"]))[0])
    )
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=arg_dict["transparent"])
    plt.show()

if __name__ == '__main__':

    transparent_flag = False
    mom = 805

    arg_dict = {
        "data": f"../results/root/htof_proton/htof_proton_eta_lambda_{mom}_hit_position.root",
        "title": r"$K^-p\rightarrow \Lambda\eta$",
        "transparent" : transparent_flag
    }
    plot(arg_dict)

    arg_dict = {
        "data": f"../results/root/htof_proton/htof_proton_kp_{mom}_hit_position.root",
        "title": r"$K^-p\rightarrow K^-p$",
        "transparent" : transparent_flag
    }
    plot(arg_dict)

    arg_dict = {
        "data": f"../results/root/htof_proton/htof_proton_pimsigmap_{mom}_hit_position.root",
        "title": r"$K^-p\rightarrow \pi^-\Sigma^+$",
        "transparent" : transparent_flag
    }
    plot(arg_dict)

    arg_dict = {
        "data": f"../results/root/htof_proton/htof_proton_pi0sigma0_{mom}_hit_position.root",
        "title": r"$K^-p\rightarrow \pi^0\Sigma^0$",
        "transparent" : transparent_flag
    }
    plot(arg_dict)

    arg_dict = {
        "data": f"../results/root/htof_proton/htof_proton_pi0lambda_{mom}_hit_position.root",
        "title": r"$K^-p\rightarrow \pi^0\Lambda$",
        "transparent" : transparent_flag
    }
    plot(arg_dict)

