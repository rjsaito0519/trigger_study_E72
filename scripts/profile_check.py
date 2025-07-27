import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys
import lmfit.models as lfm

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

def gauss_fit(values, edges_x):
    width = edges_x[1] - edges_x[0]
    x = (edges_x + np.full_like(edges_x, width/2.))[:-1]
    data = np.sum(values.T, axis=0)

    model = lfm.GaussianModel()
    params = model.guess(x = x, data = data)
    result = model.fit(x = x, data = data, params = params)
    # print(result.fit_report())
    # result.plot()
    # plt.show()

    return result.params["center"].value 

def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../data/tmp/{}".format(arg_dict["data"]))
    file = uproot.open(root_file_path)


    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    # -- hit pattern -----
    values, edges_x, edges_y = file["BAC_profile"].to_numpy()  # ヒストグラムのデータを取得
    masked_values = np.ma.masked_where(values == 0, values)       # 値が0またはNaNの部分をマスクする
    mesh = ax1.pcolormesh(
        edges_x, edges_y, masked_values.T,
        shading="flat", cmap='viridis'
    )
    center = gauss_fit(values, edges_x)
    ax1.axvline(center, ls = "dashed", color = "red")

    x_bac = 115.
    y_bac = 115.
    offset = arg_dict["offset"]
    rec = (-x_bac/2.+offset, -y_bac/2.)
    ax1.add_patch( patches.Rectangle(xy=rec, width=x_bac, height=y_bac, ec='C1', lw = 3, fill=False) )
    ax1.text(0.7, 0.5, "BAC\n({:.1f} mm)".format(center), 
        transform=ax1.transAxes,  # 軸座標を基準にする
        ha="left", va="center", color="k", zorder=2, fontsize = 30) 




    values, edges_x, edges_y = file["TGT_profile"].to_numpy()  # ヒストグラムのデータを取得
    masked_values = np.ma.masked_where(values == 0, values)       # 値が0またはNaNの部分をマスクする
    mesh = ax2.pcolormesh(
        edges_x, edges_y, masked_values.T,  
        shading="flat", cmap='viridis'
    )
    center = gauss_fit(values, edges_x)
    ax2.axvline(center, ls = "dashed", color = "red")
    ax2.text(0.7, 0.5, "Target\n({:.1f} mm)".format(center), 
        transform=ax2.transAxes,  # 軸座標を基準にする
        ha="left", va="center", color="k", zorder=2, fontsize = 30) 



    values, edges_x, edges_y = file["KVC_profile"].to_numpy()  # ヒストグラムのデータを取得
    masked_values = np.ma.masked_where(values == 0, values)       # 値が0またはNaNの部分をマスクする
    mesh = ax3.pcolormesh(
        edges_x, edges_y, masked_values.T,  
        shading="flat", cmap='viridis'
    )
    # center = gauss_fit(values, edges_x)
    # ax3.axvline(center, ls = "dashed", color = "red")
    ax3.text(0.15, 0.5, "KVC", 
        transform=ax3.transAxes,  # 軸座標を基準にする
        ha="left", va="center", color="k", zorder=2, fontsize = 30) 

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(-100, 250)
    
    ax1.set_ylim(-89, 89)
    ax1.set_xticklabels([])
    ax2.set_ylim(-60, 60)
    ax2.set_xticklabels([])
    ax3.set_ylim(-80, 80)

    fig.supxlabel("X [mm]", x = (0.15 + 0.98)/2.)
    fig.supylabel("Y [mm]")
    fig.suptitle(arg_dict["title"], x = (0.15 + 0.98)/2.)
        

    plt.subplots_adjust(left = 0.15, right=0.98, top=0.93, bottom = 0.12, wspace=0.02, hspace=0.01)


    img_save_path = os.path.join(
        script_dir, 
        "../results/img/profile_{}.png".format(os.path.splitext(os.path.basename(arg_dict["data"]))[0])
    )
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=300, transparent=True)
    plt.show()


if __name__ == '__main__':

    # -- 735 -----
    arg_dict = {
        "data": "original_pos.root",
        "offset": -14.,
        "title": r"offset = $-14$ mm (735 MeV/c)"
    }
    plot(arg_dict)

    arg_dict = {
        "data": "no_shift.root",
        "offset": 0.,
        "title": r"offset = $0.$ mm (735 MeV/c)"
    }
    plot(arg_dict)

    
    # -- 645 -----
    arg_dict = {
        "data": "original_pos_645.root",
        "offset": -14.,
        "title": r"offset = $-14$ mm (645 MeV/c)"
    }
    plot(arg_dict)

    arg_dict = {
        "data": "no_shift_645.root",
        "offset": 0.,
        "title": r"offset = $0.$ mm (645 MeV/c)"
    }
    plot(arg_dict)



    # -- 785 -----
    arg_dict = {
        "data": "original_pos_785.root",
        "offset": -14.,
        "title": r"offset = $-14$ mm (785 MeV/c)"
    }
    plot(arg_dict)

    arg_dict = {
        "data": "no_shift_785.root",
        "offset": 0.,
        "title": r"offset = $0.$ mm (785 MeV/c)"
    }
    plot(arg_dict)
    