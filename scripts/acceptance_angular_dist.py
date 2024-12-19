import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = False
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ
script_dir = os.path.dirname(os.path.abspath(__file__))


def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values*100

def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["data"]))
    file = uproot.open(root_file_path)

    center, edge, value_all = get_hist_data(file, "acceptance")
    _, _, value_mp2 = get_hist_data(file, "acceptance_mp2")
    _, _, value_htofp = get_hist_data(file, "acceptance_htofp")

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    ax.hist(center, bins=edge, weights=value_all,   lw = 2, histtype='step', color="k", label = "HTOF Mp2 || Forward Proton", zorder = 3)
    ax.hist(center, bins=edge, weights=value_mp2,   lw = 2, histtype='step', color="C0", label = "HTOF Mp2")
    ax.hist(center, bins=edge, weights=value_htofp, lw = 2, histtype='step', color="C1", label = "Forward Proton")
    ax.set_ylim(0, 100)

    ax.set_xlabel(r"$\cos \theta_{\rm CM}$")
    ax.set_ylabel("Acceptance [%]")
    ax.set_title(arg_dict["title"])
    plt.subplots_adjust(left = 0.15, right=0.98, top=0.92, bottom = 0.12)

    img_save_path = os.path.join(
        script_dir, 
        "../results/img/acceptance/{}.pdf".format(os.path.splitext(os.path.basename(arg_dict["data"]))[0])
    )
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=400, transparent=True)

    # -- for legend plot -----
    fig_legend = plt.figure(figsize=(8, 8))  # 凡例専用のFigureを作成
    fig_legend.legend(
        handles=ax.get_legend_handles_labels()[0],  # ハンドル（図の要素）
        labels=ax.get_legend_handles_labels()[1],   # ラベル（凡例テキスト）
        loc='center',  # 凡例の位置
    )    
    fig_legend.savefig(os.path.join(script_dir, "../results/img/acceptance/legend.pdf"), format='pdf', bbox_inches='tight', dpi=400, transparent=True)
    plt.close(fig_legend)

    plt.show()

if __name__ == '__main__':
    
    # +------------+
    # | eta Lambda |
    # +------------+
    arg_dict = {
        "data": "acceptance_mom735_flat_eta_lambda_2212.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
    }
    plot(arg_dict)

    
    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "data": "acceptance_mom735_flat_kp_0.root",
        "title": r"$K^-p\rightarrow K^-p$",
    }
    plot(arg_dict)
    
    # +-----+
    # | K0n |
    # +-----+
    arg_dict = {
        "data": "acceptance_mom735_flat_k0n_211.root",
        "title": r"$K^-p\rightarrow \overline{K}^0_sn$ ($\overline{K}^0_s \rightarrow \pi^+\pi^-$)",
    }
    plot(arg_dict)
    
    # +------------+
    # | pi+ Sigma- |
    # +------------+
    arg_dict = {
        "data": "acceptance_mom735_flat_pi+sigma-_2112.root",
        "title": r"$K^-p\rightarrow \pi^+\Sigma^-$ ($\Sigma^- \rightarrow n\pi^-$)",
    }
    plot(arg_dict)

    # +------------+
    # | pi- Sigma+ |
    # +------------+
    arg_dict = {
        "data": "acceptance_mom735_flat_pi-sigma+_2212.root",
        "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow p\pi^0$)",
    }
    plot(arg_dict)

    arg_dict = {
        "data": "acceptance_mom735_flat_pi-sigma+_2112.root",
        "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow n\pi^+$)",
    }
    plot(arg_dict)
    

    # +------------+
    # | pi0 Sigma0 |
    # +------------+
    arg_dict = {
        "data": "acceptance_mom735_flat_pi0sigma0_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Sigma^0$ ($\Sigma^0 \rightarrow \Lambda\gamma \rightarrow p\pi^-\gamma$)",
    }
    plot(arg_dict)


    # +------------+
    # | pi0 Lambda |
    # +------------+
    arg_dict = {
        "data": "acceptance_mom735_flat_pi0lambda_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
    }
    plot(arg_dict)
