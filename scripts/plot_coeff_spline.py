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
plt.rcParams['axes.grid'] = True
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ
script_dir = os.path.dirname(os.path.abspath(__file__))


def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../data/spline/{}".format(arg_dict["data"]))
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    ax_list = [ax1, ax2, ax3]
    for i, ax in enumerate(ax_list):
        x = tree[f"A{i}_spline_x"][0]
        y = tree[f"A{i}_spline_y"][0]        
        ax.plot(x, y)

        point_x = tree[f"A{i}_coeff_x"][0]
        point_y = tree[f"A{i}_coeff_y"][0]
        ax.plot(point_x, point_y, "o", ms = 5, color = "C3")

        if (ax != ax3):
            ax.set_xticklabels([])


    ax1.set_ylabel("$A_0$")
    ax2.set_ylabel("$A_1$")
    ax3.set_ylabel("$A_2$")
    ax3.set_xlabel(r"$p_{K}^{\rm Lab}$ [MeV/c]")
    ax1.set_title(arg_dict["title"])

    plt.subplots_adjust(left = 0.13, right=0.98, top=0.92, bottom = 0.12, hspace=0.02)
    img_save_path = os.path.join(script_dir, "../results/img/yield/example_{}.pdf".format(os.path.splitext( arg_dict["data"] )[0]))
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

    plt.show()


if __name__ == '__main__':
    
    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "data": "Kp_bubble1970_spline.root",
        "title": r"$K^-p\rightarrow K^-p$",
    }
    plot(arg_dict)