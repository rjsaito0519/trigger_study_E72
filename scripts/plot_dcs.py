# Standard imports
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import uproot

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

mass_Lambda  = 1115.683
mass_eta     =  547.862
mass_kaon    =  493.677
mass_proton  =  938.27208816
mass_neutron = 939.56542052

mass_tot = mass_Lambda + mass_eta
range_left  = mass_tot - 75
range_right = mass_tot + 75


def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["rootdata"]))
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    n_x = arg_dict["n_x"]
    n_y = arg_dict["n_y"]

    fig, axes = plt.subplots(n_y, n_x, figsize=(10, 12))

    ylim_list = []
    tmp_ylim = []

    for i, ax in enumerate(axes.flat):
        if i >= len(tree["mom"]):
            ax.axis("off")
        else:
            ax.errorbar(
                tree["cos_theta"][i], tree["diff_cs"][i], 
                xerr = tree["cos_theta_err"][i], 
                yerr = tree["diff_cs_err"][i], 
                fmt = "s", 
                capsize = 0, 
                markeredgecolor = "red", 
                ms = 0, 
                ecolor='red', 
                color="red", 
                markeredgewidth = 1.5, 
                zorder = 3
            )
            ax.grid(color="gray", linestyle="dotted", linewidth=1)
            ax.set_xlim(-1.2, 1.2)
            ax.text(
                0.5, 0.85, "{:.1f} MeV/c".format(tree["mom"][i]), 
                transform=ax.transAxes,
                ha='center', va='center',
                fontsize = 12
            )

            if i < len(tree["mom"])-n_x:
                ax.set_xticklabels([])
            if i%n_x!=0:
                ax.set_yticklabels([])

            if i != 0 and i%n_x==0:
                ylim_list.append(np.max(tmp_ylim))
                tmp_ylim = []
            else:
                tmp_ylim.append(np.max(tree["diff_cs"][i]))
    ylim_list.append(np.max(tmp_ylim))

    index = 0
    for i, ax in enumerate(axes.flat):
        if i >= len(tree["mom"]):
            break
        
        offset = ylim_list[index]*0.1
        ax.set_ylim(-offset, ylim_list[index]+offset)
        if i != 0 and i%n_x==0:
            index += 1
    
    plt.subplots_adjust(left = 0.12, right=0.99, top=0.93, bottom = 0.1, hspace=0.01, wspace=0.01)
    fig.suptitle(arg_dict["title"], x = 0.12+(0.99-0.12)/2)
    fig.supxlabel(r"$\cos \theta_{\rm CM}$", x = 0.12+(0.99-0.12)/2)
    fig.supylabel(r"d$\sigma$/d$\Omega$ [mb/sr]")
    
    img_path = os.path.join(script_dir, "../results/img/yield/{}".format(arg_dict["savename"]))
    plt.savefig(img_path, format='pdf', dpi=600, transparent=True)
    plt.show()


if __name__ == '__main__':

    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "rootdata" : "dcs_kp_9999.root",
        "n_x" : 6,
        "n_y" : 7,
        "title": r"$K^-p\rightarrow K^-p$ (E72 Expected result)",
        "savename": "kp_legendre_result.pdf"
    }
    plot(arg_dict)

    # +-----+
    # | K0n |
    # +-----+
    arg_dict = {
        "rootdata" : "dcs_k0n_211.root",
        "n_x" : 6,
        "n_y" : 10,
        "title": r"$K^-p \rightarrow \overline{K}^0n$ (E72 Expected result)",
        "savename": "k0n_legendre_result.pdf"
    }
    plot(arg_dict)

    # +---------------+
    # | Kp->pi+Sigma- |
    # +---------------+
    arg_dict = {
        "rootdata" : "dcs_pi+sigma-_2112.root",
        "n_x" : 6,
        "n_y" : 10,
        "title": r"$K^-p \rightarrow \pi^+\Sigma^-$ (E72 Expected result)",
        "savename": "pi+sigma-_legendre_result.pdf"    
    }
    plot(arg_dict)