import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys
import lmfit as lf
import lmfit.models as lfm
import pprint

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

# bac_optimize_mom645_-12.5.root   

data = {
    645: [],
    735: [],
    933: []
}
for mom in [645, 735, 933]:
    for pos in np.arange(-5, -25.1, -0.5):
        root_file_path = os.path.join(script_dir, f"../results/root/bac_pos_optimize/bac_optimize_mom{mom:.0f}_{pos:.1f}.root")
        file = uproot.open(root_file_path)
        tree = file["tree"].arrays(library="np")
        ratio = tree["n_hit_kaon"][0]/tree["n_kaon"][0]
        data[mom].append([pos, ratio, tree["n_hit_tgt"][0]])
    data[mom] = np.array(data[mom])

pprint.pprint(data)


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

for i, mom in enumerate([645, 735, 933]):
    x = data[mom][:, 0]
    y = data[mom][:, 2]/1000000
    ax.plot(x, y, "--o", color = f"C{i}", label = r"$p_K = ${} [MeV/c]".format(mom))
    
    model = lfm.GaussianModel()
    params = model.guess(x = x, data = y)

    result = model.fit(x = x, data = y, params=params, method='leastsq')
    # print(result.fit_report())
    fit_x = np.linspace(-0.5, -25.5, 10000)
    fit_y = result.eval_components(x=fit_x)["gaussian"]
    ax.plot(fit_x, fit_y, "--", color = f"C{i}")

    # ax.axvline(result.params["center"].value, ls = "dashed", color = f"C{i}")

    print(result.params["center"].value)

ax.set_ylabel(r"$N_{\text{hit TGT}}\ /\ N_{\text{generate}}$")
ax.set_xlabel("Baem and BAC offset [mm]")

plt.legend(fontsize = 20)
plt.subplots_adjust(left = 0.17, right=0.98, top=0.98, bottom = 0.12)

# img_save_path = os.path.join( script_dir, "../results/img/cvc.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)

plt.show()