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

root_file_path = os.path.join(script_dir, "../results/root/kvc_profile.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)

i = 0
for mom in [645, 735, 805]:
    idx = tree["mom"] == mom
    ax.errorbar(
        tree["z_pos"][idx],
        tree["mean_val"][idx],
        yerr = tree["mean_err"][idx],
        fmt = "o", capsize = 0, markeredgecolor = "k", ms = 8, 
        ecolor='k', color=f'C{i}', markeredgewidth = 0.1, zorder = 3,
        label = r"$p_K = ${} MeV/c".format(mom)
    )
    model = lfm.LinearModel()
    params = model.guess(x = tree["z_pos"][idx], data = tree["mean_val"][idx])
    result = model.fit(x = tree["z_pos"][idx], data = tree["mean_val"][idx], weights = 1.0/tree["mean_err"][idx], params=params, method='leastsq')
    print(mom)
    print(result.fit_report())
    print("------------------")
    fit_x = np.linspace(min(tree["z_pos"][idx]), max(tree["z_pos"][idx]), 10)
    fit_y = result.eval_components(x=fit_x)["linear"]
    ax.plot(fit_x, fit_y, color = f"C{i}", ls = "dashed")
    i += 1


plt.legend(fontsize = 18)
ax.set_xlabel("Z [mm]")
ax.set_ylabel("Beam Center [mm]")
plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.12)
plt.show()

# ax1.yaxis.set_major_formatter(ptick.EngFormatter())
# ax1.set_xlabel(r"thickness [mm]")


# # img_save_path = os.path.join(script_dir, "../results/img/yield/beam_distribution_scan.pdf")
# # os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# # plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

# plt.show()