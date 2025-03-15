import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import numpy as np
import uproot
import os
import sys
import matplotlib.gridspec as gridspec


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


root_file_path = os.path.join(script_dir, "../results/root/kvc_pos_optimize.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

beam_x = tree["beam_x"][0]/10
beam_y = tree["beam_y"][0]/10
beam_z = tree["beam_z"][0]/10
beam_u = tree["beam_u"][0]
beam_v = tree["beam_v"][0]

proton_x = tree["proton_x"][0]/10
proton_y = tree["proton_y"][0]/10
proton_z = tree["proton_z"][0]/10
proton_u = tree["proton_u"][0]
proton_v = tree["proton_v"][0]

fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
ax4 = fig.add_subplot(gs[1, 2])
ax5 = fig.add_subplot(gs[1, 3])

beam_profile = [
    [],
    [],
    [],
    [],
]
proton_profile = [
    [],
    [],
    [],
    [],
]

z = np.array([400, 450, 550, 650, 750, 800])/10
for i in range(len(beam_x)):
    x = beam_x[i] + beam_u[i]*(z-beam_z[i])
    ax1.plot(z, x, color = "C0", alpha = 0.2, lw = 0.1)
    for j in range(4):
        beam_profile[j].append(x[j+1])

for i in range(len(proton_x)):
    x = proton_x[i] + proton_u[i]*(z-proton_z[i])
    ax1.plot(z, x, color = "C1", alpha = 0.2, lw = 0.1)
    for j in range(4):
        proton_profile[j].append(x[j+1])

y_min, y_max = np.inf, -np.inf
mpv = []
for i, ax in enumerate([ax2, ax3, ax4, ax5]):
    hist_data = ax.hist(beam_profile[i], bins = np.linspace(-70, 50, 101), histtype='step', color="C0")
    index = np.argmax(hist_data[0])
    mpv.append(hist_data[1][index])
    ax.hist(proton_profile[i], bins = np.linspace(-70, 50, 101), histtype='step', color="C1")
    ax.set_yticklabels([])
    tmp_y_min, tmp_y_max = ax.get_ylim()
    if tmp_y_min < y_min:
        y_min = tmp_y_min
    if tmp_y_max > y_max:
        y_max = tmp_y_max
            
    ax1.axvline(45+10*i, ls = "dashed", color = "red")

ax2.set_ylim(y_min, y_max)
ax3.set_ylim(y_min, y_max)
ax4.set_ylim(y_min, y_max)
ax5.set_ylim(y_min, y_max)

ax2.fill_betweenx([y_min, y_max], 16-2.6*4, 16+2.6*4, fc="C2", alpha = 0.2, ec = "None", zorder = 0)
ax3.fill_betweenx([y_min, y_max], 16-2.6*4+(mpv[1]-mpv[0]), 16+2.6*4+(mpv[1]-mpv[0]), fc="C2", alpha = 0.2, ec = "None", zorder = 0)
ax4.fill_betweenx([y_min, y_max], 16-2.6*4+(mpv[2]-mpv[0]), 16+2.6*4+(mpv[2]-mpv[0]), fc="C2", alpha = 0.2, ec = "None", zorder = 0)
ax5.fill_betweenx([y_min, y_max], 16-2.6*4+(mpv[3]-mpv[0]), 16+2.6*4+(mpv[3]-mpv[0]), fc="C2", alpha = 0.2, ec = "None", zorder = 0)

for i, ax in enumerate([ax2, ax3, ax4, ax5]):
    ax.text(0.1, 0.9, r"$Z = {}$ cm".format(45+10.0*i) + "\n" + r"$X = {:.1f}$ cm".format(16+(mpv[i]-mpv[0])), 
        transform=ax.transAxes,  # 軸座標を基準にする
        ha="left", va="center", color="k", zorder=2, fontsize = 20) 

ax1.set_ylabel("X [cm]")
ax1.set_xlabel("Z [cm]")
fig.supxlabel("X [cm]", x = (0.12+0.98)/2.)
plt.subplots_adjust(left = 0.12, right=0.98, top=0.98, bottom = 0.12, wspace = 0.02, hspace = 0.25)
img_save_path = os.path.join(script_dir, "../results/img/track.jpg")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='jpg', bbox_inches='tight', dpi=600, transparent=True)

plt.show()

# # +-----------------------------+
# # | scan Kaon beam distribution |
# # +-----------------------------+
# fig = plt.figure(figsize=(10, 8))
# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)

# # ax1.set_xlim(-5, 410)
# # ax2.set_xlim(-5, 410)

# ax1.set_xlim(-2, 50)
# ax2.set_xlim(-2, 50)

# # label_list1 = ["kaon_thin", "electron_thin", "kaon_thick", "electron_thick"]
# label_list1 = ["kaon_thin_woair", "electron_thin_noair", "kaon_thick_noair", "electron_thick_noair"]
# color_list1 = [
#     adjust_lightness(mcolors.to_rgb('C0'), 1.),
#     adjust_lightness(mcolors.to_rgb('C1'), 1.),
#     adjust_lightness(mcolors.to_rgb('C0'), 0.8),
#     adjust_lightness(mcolors.to_rgb('C1'), 0.8),
# ]

# # label_list2 = ["kaon_1seg_thin", "electron_1seg_thin", "kaon_1seg_thick", "electron_1seg_thick"]
# label_list2 = ["kaon_1seg_thin_woair", "electron_1seg_thin_noair", "kaon_1seg_thick", "electron_1seg_thick_noair"]
# color_list2 = [
#     adjust_lightness(mcolors.to_rgb('C0'), 1.),
#     adjust_lightness(mcolors.to_rgb('C1'), 1.),
#     adjust_lightness(mcolors.to_rgb('C0'), 0.8),
#     adjust_lightness(mcolors.to_rgb('C1'), 0.8),
# ]

# for label, color in zip(label_list1, color_list1):
#     root_file_path = os.path.join(script_dir, f"../results/root/{label}.root")
#     file = uproot.open(root_file_path)
#     center, edge, value = get_hist_data(file, f"hit")
#     ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color=color, zorder = 3, label = label)

# for label, color in zip(label_list2, color_list2):
#     root_file_path = os.path.join(script_dir, f"../results/root/{label}.root")
#     file = uproot.open(root_file_path)
#     center, edge, value = get_hist_data(file, f"hit")
#     ax2.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color=color, zorder = 3, label = label)


# # ax1.yaxis.set_major_formatter(ptick.EngFormatter())

# ax1.set_xticklabels([])

# ax2.set_xlabel(r"$N_{\rm p.e.}$")
# # ax1.set_ylabel("")
# plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.12, hspace = 0)

# ax1.legend(fontsize = 20)
# ax2.legend(fontsize = 20)


# img_save_path = os.path.join(script_dir, "../results/img/npe.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)

# plt.show()
