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


root_file_path = os.path.join(script_dir, "../data/n_kaon.root")
file = uproot.open(root_file_path)

# +-----------------------------+
# | scan Kaon beam distribution |
# +-----------------------------+
fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(111)
sum_n_kaon_scan = 0

for i, mom in enumerate([685, 705, 725, 745, 765]):
    center, edge, value = get_hist_data(file, f"mom{mom}")
    sum_n_kaon_scan += np.sum(value)
    width = (edge[1] - edge[0])
    ax1.hist(center, bins=edge, weights=value/width, lw = 1.5, histtype='step', color=f"C{i}", zorder = 3, label = f"{mom:.0f} MeV/c")

center, edge, value = get_hist_data(file, "mom_scan")
ax1.hist(center, bins=edge, weights=value/width, lw = 1.5, histtype='step', color="gray", zorder = 3)
ax1.text(0.07, 0.93, f"Total = {sum_n_kaon_scan/10**9:.1f}" + r" $\times 10^9$ Kaons", 
        transform=ax1.transAxes,  # 軸の座標系を基準
        fontsize=28, 
        color='k')

ax1.yaxis.set_major_formatter(ptick.EngFormatter())
ax1.set_xlabel(r"$p_{K}^{\rm Lab}$ [MeV/c]")
ax1.set_ylabel("Counts [/(MeV/c)]")
ax1.set_xlim(649, 791)

plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.12)
img_save_path = os.path.join(script_dir, "../results/img/yield/beam_distribution_scan.pdf")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

plt.show()

# +------------------------+
# | Kaon beam distribution |
# +------------------------+
fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(111)

center, edge, value = get_hist_data(file, "mom735")
ax1.hist(center, bins=edge, weights=value/width, lw = 1.5, histtype='step', color="C0", zorder = 3)
ax1.text(0.07, 0.93, f"Total = {np.sum(value)/10**9:.1f}" + r" $\times 10^9$ Kaons", 
        transform=ax1.transAxes,  # 軸の座標系を基準
        fontsize=28, 
        color='k')

ax1.yaxis.set_major_formatter(ptick.EngFormatter())
ax1.set_xlabel(r"$p_{K}^{\rm Lab}$ [MeV/c]")
ax1.set_ylabel("Counts [/(MeV/c)]")
ax1.set_xlim(649, 791)

plt.subplots_adjust(left = 0.16, right=0.98, top=0.98, bottom = 0.12)
img_save_path = os.path.join(script_dir, "../results/img/yield/beam_distribution_735.pdf")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

plt.show()

print( format(sum_n_kaon_scan + np.sum(value), "E") )