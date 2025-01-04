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

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values


root_file_path = os.path.join(script_dir, "../data/angle.root")
file = uproot.open(root_file_path)

fig = plt.figure(figsize=(8, 5))
ax1 = fig.add_subplot(111)
sum_n_kaon_scan = 0

center, edge, value = get_hist_data(file, f"angle")
ax1.hist(center, bins=edge, weights=value, lw = 1., histtype='step', color="k", zorder = 3)

ax1.yaxis.set_major_formatter(ptick.EngFormatter())
ax1.set_xlabel(r"Angle [deg]")
ax1.set_xlim(12, 29)

plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.16)
img_save_path = os.path.join(script_dir, "../results/img/KVC_angle.pdf")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

plt.show()