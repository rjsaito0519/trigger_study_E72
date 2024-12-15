import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import numpy as np
import uproot
import os

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

    # fig = plt.figure(figsize=(12, 8))
    # ax  = fig.add_subplot(111)

    # # -- histogram -----
    # center, edge, value = get_hist_data(file, "hitpat_all")
    # ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 2)
    # ax.xaxis.set_major_formatter(ptick.EngFormatter())
    # ax.yaxis.set_major_formatter(ptick.EngFormatter())
    

    edep_list = [
        # hist_name       color  label
        ["edep_proton",   "C3", "proton"],
        ["edep_piplus",   "C1", r"$\pi^+$"],
        ["edep_piminus",  "C0", r"$\pi^-$"],
        ["edep_kminus",   "C2", r"$K^-$"],
        ["edep_electron", "C5", r"$e^-, e^+$"],
        ["edep_muon",     "C6", r"$\mu^-, \^mu^+$"],
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
    # ax1.legend(fontsize = 20, handletextpad = 0.5, handlelength=0.7, loc = "upper left")
    ax1.legend(fontsize = 20, handletextpad = 0.5, handlelength=0.7)

    # -- hit pattern -----
    ax2 = ax[1]
    values, edges_x, edges_y = file["hitpat_vs_edep"].to_numpy()  # ヒストグラムのデータを取得
    masked_values = np.ma.masked_where(values == 0, values)       # 値が0またはNaNの部分をマスクする
    mesh = ax2.pcolormesh(
        edges_y, edges_x, masked_values,  # xとyを逆にしているので転置は不要
        shading="flat", cmap='viridis'
    )


    # if arg_dict["log"]:
    #     counts, xedges, yedges, cbar = ax2.hist2d(edep_vs_seg[:, 0], edep_vs_seg[:, 1], bins=[np.linspace(0., 12, 121), np.linspace(-0.5, 33.5, 35)], cmap=cmap, norm=matplotlib.colors.LogNorm())
    # else:
    #     counts, xedges, yedges, cbar = ax2.hist2d(edep_vs_seg[:, 0], edep_vs_seg[:, 1], bins=[np.linspace(0., 12, 121), np.linspace(-0.5, 33.5, 35)], cmap=cmap)
    # cbar.set_clim(1, np.max(counts))
    # ax2.add_patch( patches.Rectangle(xy=rec, width=13, height=n_seg, ec='C3', lw = 3, fill=False) )
    ax2.set_ylabel("HTOF segment")
    ax2.set_xlabel("Energy deposit [MeV]")

    # #余白の調整
    plt.subplots_adjust(left = 0.12, right=0.89, top=0.91)
    plt.subplots_adjust(wspace=0.02, hspace=0.01)
    cax = plt.axes([0.895, 0.1, 0.025, 0.484])
    fig.colorbar(mesh, cax=cax, pad=0.01, fraction=0.1)
    # plt.savefig(f"./img/{name}.png", dpi=600, transparent=True)
    plt.show()



if __name__ == '__main__':
    arg_dict = {
        "data": "htof_mom735_eta_lambda.root",
        "title": r"$K^-p\rightarrow \eta\Lambda$",
        "log": False
    }
    plot(arg_dict)



# log_flag = False

# # -- eta lambda -----
# name = "mom735_eta_lambda_htof"
# title = r"$K^-p\rightarrow \eta\Lambda$"
# rec = (4.0, 17-0.5)
# n_seg = 5

# # -- Kp -----
# name = "mom735_kp_htof"
# title = r"$K^-p\rightarrow K^-p$"
# rec = (4.0, 15-0.5)
# n_seg = 11

# # -- pi- Sigma+ -----
# name = "mom735_pi-sigma+_htof"
# title = r"$K^-p\rightarrow \pi^-\Sigma^+$"
# rec = (3.0, 15-0.5)
# n_seg = 11

# # -- beam -----
# name = "mom735_beam_htof"
# title = r"$K^-$ beam"
# rec = (3.0, 6-0.5)
# n_seg = 28
# log_flag = True


# def load_data(path):
#     data = np.genfromtxt(path, skip_header=1, delimiter=",")
#     if len(data) == 0:
#         data = np.array([[np.nan, np.nan]])
#     elif len(data) == 2:
#         data = np.vstack((data, [np.nan, np.nan]))
#     return data

# data_proton = load_data(f"../.csv_data/{name}_proton.csv")
# data_piplus = load_data(f"../.csv_data/{name}_piplus.csv")
# data_piminus = load_data(f"../.csv_data/{name}_piminus.csv")
# data_kaon = load_data(f"../.csv_data/{name}_kaon.csv")
# data_electron = load_data(f"../.csv_data/{name}_electron.csv")
# data_muon = load_data(f"../.csv_data/{name}_muon.csv")

# edep = {
#     2212 : data_proton[:, 0],
#      211 : data_piplus[:, 0],
#     -211 : data_piminus[:, 0],
#     -321 : data_kaon[:, 0],
#       11 : data_electron[:, 0],
#       13 : data_muon[:, 0]
# }

# edep_vs_seg = np.vstack((
#     data_proton,
#     data_piplus,
#     data_piminus,
#     data_kaon,
#     data_electron,
#     data_muon
# ))

# fig, ax = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [2, 3]})
# ax1 = ax[0]
# ax1.hist(edep[2212], bins = np.linspace(0., 12, 121), histtype='step', color = "C3", label = "proton")
# ax1.hist(edep[211],  bins = np.linspace(0., 12, 121), histtype='step', color = "C1", label = r"$\pi^+$")
# ax1.hist(edep[-211], bins = np.linspace(0., 12, 121), histtype='step', color = "C0", label = r"$\pi^-$")
# ax1.hist(edep[-321], bins = np.linspace(0., 12, 121), histtype='step', color = "C2", label = r"$K^-$")
# ax1.hist(edep[11], bins = np.linspace(0., 12, 121), histtype='step', color = "C5", label = r"$e^-, e^+$")
# ax1.hist(edep[13], bins = np.linspace(0., 12, 121), histtype='step', color = "C6", label = r"$\mu^-, \mu^+$")

# ax1.yaxis.set_major_formatter(ptick.EngFormatter())
# ax1.set_xticks([0, 2, 4, 6, 8, 10, 12])
# ax1.axes.xaxis.set_ticklabels([])
# ax1.grid(True)
# ax1.set_axisbelow(True)
# ax1.set_title(title)
# # ax1.legend(fontsize = 20, handletextpad = 0.5, handlelength=0.7, loc = "upper left")
# ax1.legend(fontsize = 20, handletextpad = 0.5, handlelength=0.7)

# ax2 = ax[1]
# cmap = copy.copy(plt.cm.viridis)
# cmap.set_under('w', 1) # 下限以下の色を設定
# if log_flag:
#     counts, xedges, yedges, cbar = ax2.hist2d(edep_vs_seg[:, 0], edep_vs_seg[:, 1], bins=[np.linspace(0., 12, 121), np.linspace(-0.5, 33.5, 35)], cmap=cmap, norm=matplotlib.colors.LogNorm())
# else:
#     counts, xedges, yedges, cbar = ax2.hist2d(edep_vs_seg[:, 0], edep_vs_seg[:, 1], bins=[np.linspace(0., 12, 121), np.linspace(-0.5, 33.5, 35)], cmap=cmap)
# cbar.set_clim(1, np.max(counts))
# ax2.add_patch( patches.Rectangle(xy=rec, width=13, height=n_seg, ec='C3', lw = 3, fill=False) )
# ax2.set_ylabel("HTOF segment")
# ax2.set_xlabel("Energy deposit [MeV]")

# #余白の調整
# plt.subplots_adjust(left = 0.12, right=0.89, top=0.91)
# plt.subplots_adjust(wspace=0.02, hspace=0.01)
# cax = plt.axes([0.895, 0.1, 0.025, 0.484])
# fig.colorbar(cbar, cax=cax, pad=0.01, fraction=0.1)
# plt.savefig(f"./img/{name}.png", dpi=600, transparent=True)
# plt.show()