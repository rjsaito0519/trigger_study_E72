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


    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # -- meson -----
    values, edges_x, edges_y = file["meson"].to_numpy()
    masked_values = np.ma.masked_where(values == 0, values)
    mesh1 = ax1.pcolormesh(
        edges_x, edges_y/1000., masked_values.T,
        shading="flat", cmap='viridis', zorder=1
    )
    cax1 = plt.axes([0.855, 0.525, 0.03, 0.396])
    fig.colorbar(mesh1, cax=cax1, pad=0.01, fraction=0.1)
    ax1.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax1.set_xticklabels([])
    ax1.set_ylim(-0.02, 1.1)
    ax1.text(0.8, 0.9, r"scattered $K^-$", 
        transform=ax1.transAxes,  # 軸座標を基準にする
        ha="center", va="center", color="k", zorder=2)

    # -- baryon -----
    values, edges_x, edges_y = file["baryon"].to_numpy()
    masked_values = np.ma.masked_where(values == 0, values)
    mesh2 = ax2.pcolormesh(
        edges_x, edges_y/1000., masked_values.T,
        shading="flat", cmap='viridis', zorder=1
    )
    cax2 = plt.axes([0.855, 0.12, 0.03, 0.396])
    fig.colorbar(mesh2, cax=cax2, pad=0.01, fraction=0.1)
    ax2.set_ylim(-0.02, 1.1)
    ax2.text(0.2, 0.9, r"scattered $p$", 
        transform=ax2.transAxes,  # 軸座標を基準にする
        ha="center", va="center", color="k", zorder=2)

    ax1.set_ylabel(r"$p_{K}$ [GeV/c]")
    ax1.set_title(arg_dict["title"])

    ax2.set_ylabel(r"$p_{p}$ [GeV/c]")
    ax2.set_xlabel(r"$\cos \theta_{\rm CM}$")
    
    plt.subplots_adjust(left = 0.16, right=0.85, top=0.92, bottom = 0.12, hspace=0.02)

    img_save_path = os.path.join(
        script_dir, 
        "../results/img/explain/{}.png".format(os.path.splitext(os.path.basename(arg_dict["data"]))[0])
    )
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=400, transparent=True)

    plt.show()

if __name__ == '__main__':
    
    # # +------------+
    # # | eta Lambda |
    # # +------------+
    # arg_dict = {
    #     "data": "acceptance_mom735_eta_lambda_2212.root",
    #     "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
    # }
    # plot(arg_dict)

    
    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "data": "explain_mom735_flat_kp.root",
        # "data" : "explain_mom735_flat_pi-sigma+.root",
        "title": r"$K^-p\rightarrow K^-p$",
    }
    plot(arg_dict)
    
    # # +-----+
    # # | K0n |
    # # +-----+
    # arg_dict = {
    #     "data": "acceptance_mom735_k0n_211.root",
    #     "title": r"$K^-p\rightarrow \overline{K}^0_sn$ ($\overline{K}^0_s \rightarrow \pi^+\pi^-$)",
    # }
    # plot(arg_dict)
    
    # # +------------+
    # # | pi+ Sigma- |
    # # +------------+
    # arg_dict = {
    #     "data": "acceptance_mom735_pi+sigma-_2112.root",
    #     "title": r"$K^-p\rightarrow \pi^+\Sigma^-$ ($\Sigma^- \rightarrow n\pi^-$)",
    # }
    # plot(arg_dict)

    # # +------------+
    # # | pi- Sigma+ |
    # # +------------+
    # arg_dict = {
    #     "data": "acceptance_mom735_pi-sigma+_2212.root",
    #     "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow p\pi^0$)",
    # }
    # plot(arg_dict)

    # arg_dict = {
    #     "data": "acceptance_mom735_pi-sigma+_2112.root",
    #     "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow n\pi^+$)",
    # }
    # plot(arg_dict)
    

    # # +------------+
    # # | pi0 Sigma0 |
    # # +------------+
    # arg_dict = {
    #     "data": "acceptance_mom735_pi0sigma0_2212.root",
    #     "title": r"$K^-p\rightarrow \pi^0\Sigma^0$ ($\Sigma^0 \rightarrow \Lambda\gamma \rightarrow p\pi^-\gamma$)",
    # }
    # plot(arg_dict)


    # # +------------+
    # # | pi0 Lambda |
    # # +------------+
    # arg_dict = {
    #     "data": "acceptance_mom735_pi0lambda_2212.root",
    #     "title": r"$K^-p\rightarrow \pi^0\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
    # }
    # plot(arg_dict)
