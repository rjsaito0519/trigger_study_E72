import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import numpy as np
import uproot
import os
import lmfit as lf
import lmfit.models as lfm
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

def x_pos(mom, z):
    a = {
        645: 0.45015376,
        735: 0.38843815,
        805: 0.35118881
    }
    b = {
        645: -10.4413093,
        735: -8.74197715,
        805: -7.74282775
    }

    return a[mom]*z + b[mom]

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

root_file_path = os.path.join(script_dir, "../results/root/kvc_profile.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

i = 0
for mom, ax in zip([645, 735, 805], [ax1, ax2, ax3]):
    for i, z_pos in enumerate(range(450, 850, 50)):
        center, edge, value = get_hist_data(file, f"mom{mom}_{z_pos}")
        ax.hist(center-x_pos(mom, z_pos), bins=edge-x_pos(mom, z_pos), weights=value, lw = 1.5, histtype='step', color=f"C{i}", zorder = 3, label = f"{z_pos:.0f} [mm]")
    ax.set_xlim(-149, 149)
    ax.axvline(-26*4, ls = "dashed", color = "red")
    ax.axvline(+26*4, ls = "dashed", color = "red")
    ax.yaxis.set_major_formatter(ptick.EngFormatter())
    ax.text(0.05, 0.85, r"$p_K$" + f" = {mom} MeV/c",
        transform=ax.transAxes,
        fontsize=20, va="center", ha="left",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none")    
    )


for ax in [ax1, ax2]:
    ax.set_xticklabels([])

ax1.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1), title = "KVC z-pos", title_fontsize = 20)
fig.supxlabel(r"$x$ - mean [mm]")
plt.subplots_adjust(left = 0.1, right=0.78, top=0.98, bottom = 0.12, hspace=0.01)
plt.show()

# ax1.yaxis.set_major_formatter(ptick.EngFormatter())
# ax1.set_xlabel(r"thickness [mm]")


# # img_save_path = os.path.join(script_dir, "../results/img/yield/beam_distribution_scan.pdf")
# # os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# # plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

# plt.show()