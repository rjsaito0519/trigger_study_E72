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
    root_file_path1 = os.path.join(script_dir, "../data/{}".format(arg_dict["data1"]))
    file1 = uproot.open(root_file_path1)

    root_file_path2 = os.path.join(script_dir, "../data/{}".format(arg_dict["data2"]))
    file2 = uproot.open(root_file_path2)

    # +-----------------------------+
    # | scan Kaon beam distribution |
    # +-----------------------------+
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(111)

    center, edge, value = get_hist_data(file1, "mom_dist_trig")
    _, _, value_scan = get_hist_data(file2, "mom_dist_trig")
    width = (edge[1] - edge[0])
    ax1.hist(center, bins=edge, weights=value/width, lw = 1.5, histtype='step', color="k", zorder = 3, label = "All")
    ax1.hist(center, bins=edge, weights=value_scan/width, lw = 1.5, histtype='step', color="C0", zorder = 1, label = "Scan")
    ax1.hist(center, bins=edge, weights=(value-value_scan)/width, lw = 1.5, histtype='step', color="C1", zorder = 2, label = "735 MeV/c")

    ax1.yaxis.set_major_formatter(ptick.EngFormatter())
    ax1.set_xlabel(r"$p_{K}^{\rm Lab}$ [MeV/c]")
    ax1.set_ylabel("Expected Yields [/(MeV/c)]")
    # ax1.axvline(723.2931, ls = "dashed", color="red")
    ax1.set_xlim(721, 775)
    ax1.set_title(arg_dict["title"])
    ax1.legend()

    plt.subplots_adjust(left = 0.15, right=0.98, top=0.92, bottom = 0.12)
    img_save_path = os.path.join(script_dir, "../results/img/yield/yield_distribution_LambdaEta.pdf")
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

    plt.show()


if __name__ == '__main__':
    
    # +------------+
    # | eta Lambda |
    # +------------+
    arg_dict = {
        "data1": "yield_eta_lambda_2212.root",
        "data2": "yield_eta_lambda_2212_onlyscan.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$ ($\Lambda \rightarrow p\pi^-$)",
    }
    plot(arg_dict)