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
    return bin_centers, bin_edges, bin_values

def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["data"]))
    file = uproot.open(root_file_path)

    prefix = ""
    if arg_dict["mp1"]:
        prefix = "multi1_"
    
    edep_list = [
        # hist_name       color  label
        [f"{prefix}edep_proton",   "C3", "proton"],
        [f"{prefix}edep_piplus",   "C1", r"$\pi^+$"],
        [f"{prefix}edep_piminus",  "C0", r"$\pi^-$"],
        [f"{prefix}edep_kminus",   "C2", r"$K^-$"],
        [f"{prefix}edep_electron", "C5", r"$e^-, e^+$"],
        [f"{prefix}edep_muon",     "C6", r"$\mu^-, \mu^+$"],
    ]

    fig, ax = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [2, 3]})

    # -- edep -----
    ax1 = ax[0]
    for keys in edep_list:
        center, edge, value = get_hist_data(file, keys[0])
        ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color=keys[1], label = keys[2])

    ax1.yaxis.set_major_formatter(ptick.EngFormatter())
    ax1.set_xticks([0, 2, 4, 6, 8, 10, 12])
    ax1.set_xticklabels([])
    ax1.grid(True)
    ax1.set_axisbelow(True)
    ax1.set_title(arg_dict["title"])
    ax1.legend(fontsize = 20, handletextpad = 0.5, handlelength=0.7, loc = "upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

    # -- hit pattern -----
    ax2 = ax[1]
    values, edges_x, edges_y = file[f"{prefix}hitpat_vs_edep"].to_numpy()  # ヒストグラムのデータを取得
    masked_values = np.ma.masked_where(values == 0, values)       # 値が0またはNaNの部分をマスクする

    if arg_dict["log"]:
        mesh = ax2.pcolormesh(
            edges_y, edges_x, masked_values,  # xとyを逆にしているので転置は不要
            shading="flat", cmap='viridis', norm=LogNorm(vmin=1, vmax=np.max(masked_values.flatten()))
        )
    else:
        mesh = ax2.pcolormesh(
            edges_y, edges_x, masked_values,  # xとyを逆にしているので転置は不要
            shading="flat", cmap='viridis'
        )

    if arg_dict["rec"]:
        rec = (3.0, 10-0.5)
        n_seg = 20
        ax1.axvline(3.0, ls = "dashed", color = "C3", lw = 2)
        ax2.add_patch( patches.Rectangle(xy=rec, width=13, height=n_seg, ec='C3', lw = 3, fill=False) )

    ax1.set_xlim(0, 11.9)
    ax2.set_xlim(0, 11.9)
    ax2.set_ylabel("HTOF segment")
    ax2.set_xlabel("Energy deposit [MeV]")

    plt.subplots_adjust(left = 0.11, right=0.85, top=0.93, bottom = 0.1, wspace=0.02, hspace=0.01)
    cax = plt.axes([0.855, 0.1, 0.03, 0.495])
    fig.colorbar(mesh, cax=cax, pad=0.01, fraction=0.1)


    img_save_path = os.path.join(
        script_dir, 
        "../results/img/htof/{}{}.png".format(prefix, os.path.splitext(os.path.basename(arg_dict["data"]))[0])
    )
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=300, transparent=True)
    plt.show()


if __name__ == '__main__':

    arg_dict = {
        "data": "htof_mom735_beam_9999.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
        "mp1": False,
        "log": True,
        "rec": True
    }
    plot(arg_dict)
    # sys.exit()

    # +------------+
    # | eta Lambda |
    # +------------+
    arg_dict = {
        "data": "htof_mom735_eta_lambda_2212.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
        "mp1": False,
        "log": False,
        "rec": False
    }
    plot(arg_dict)

    arg_dict = {
        "data": "htof_mom735_eta_lambda_2212.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
        "mp1": True,
        "log": False,
        "rec": False
    }
    plot(arg_dict)
    

    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "data": "htof_mom735_kp_0.root",
        "title": r"$K^-p\rightarrow K^-p$",
        "mp1": False,
        "log": False,
        "rec": True
    }
    plot(arg_dict)

    arg_dict = {
        "data": "htof_mom735_kp_0.root",
        "title": r"$K^-p\rightarrow K^-p$",
        "mp1": True,
        "log": False,
        "rec": True
    }
    plot(arg_dict)    
    

    # +------------+
    # | pi- Sigma+ |
    # +------------+
    arg_dict = {
        "data": "htof_mom735_pi-sigma+_2212.root",
        "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow p\pi^0$)",
        "mp1": False,
        "log": False,
        "rec": True
    }
    plot(arg_dict)

    arg_dict = {
        "data": "htof_mom735_pi-sigma+_2212.root",
        "title": r"$K^-p\rightarrow \pi^-\Sigma^+$ ($\Sigma^+ \rightarrow p\pi^0$)",
        "mp1": True,
        "log": False,
        "rec": True
    }
    plot(arg_dict)


    # +------------+
    # | pi0 Sigma0 |
    # +------------+
    arg_dict = {
        "data": "htof_mom735_pi0sigma0_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Sigma^0$ ($\Sigma^0 \rightarrow \Lambda\gamma \rightarrow p\pi^-\gamma$)",
        "mp1": False,
        "log": False,
        "rec": True
    }
    plot(arg_dict)

    arg_dict = {
        "data": "htof_mom735_pi0sigma0_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Sigma^0$ ($\Sigma^0 \rightarrow \Lambda\gamma \rightarrow p\pi^-\gamma$)",
        "mp1": True,
        "log": False,
        "rec": True
    }
    plot(arg_dict)


    # +------------+
    # | pi0 Lambda |
    # +------------+
    arg_dict = {
        "data": "htof_mom735_pi0lambda_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
        "mp1": False,
        "log": False,
        "rec": True
    }
    plot(arg_dict)

    arg_dict = {
        "data": "htof_mom735_pi0lambda_2212.root",
        "title": r"$K^-p\rightarrow \pi^0\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
        "mp1": True,
        "log": False,
        "rec": True
    }
    plot(arg_dict)
