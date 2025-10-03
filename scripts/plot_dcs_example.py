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

def legendre(x, a0, a1, a2, a3, a4, a5):
    term0 = 1
    term1 = x
    term2 = 3/2 * x**2 - 1/2
    term3 = 5/2 * x**3 - 3/2 * x
    term4 = 35/8 * x**4 - 30/8 * x**2 + 3/8
    term5 = 63/8 * x**5 - 70/8 * x**3 + 15/8 * x
    return a0*term0 + a1*term1 + a2*term2 + a3*term3 + a4*term4 + a5*term5


def plot(arg_dict):
    root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["rootdata"]))
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111)

    i = arg_dict["idx"]

    ax.errorbar(
        tree["cos_theta"][i], tree["diff_cs"][i], 
        xerr = tree["cos_theta_err"][i], 
        yerr = tree["diff_cs_err"][i], 
        fmt = "o", 
        capsize = 0, 
        markeredgecolor = "red", 
        ms = 5, 
        ecolor='red', 
        color="w", 
        markeredgewidth = 1.5, 
        zorder = 3,
    )
    ax.grid(color="gray", linestyle="dotted", linewidth=1)
    ax.set_xlim(-1.2, 1.2)

    ax.text(
        0.5, 0.9, "{:.1f} MeV/c".format(tree["mom"][i]), 
        transform=ax.transAxes,
        ha='center', va='center',
        fontsize = 40,
        bbox=dict(facecolor="yellow", alpha=0.5, edgecolor="none")
    )

    fit_x = np.linspace(-1.0, 1.0, 1000)
    if arg_dict["coeff_flag"]:
        fit_y = legendre(fit_x, tree["coeff"][i][0], tree["coeff"][i][1], tree["coeff"][i][2], tree["coeff"][i][3], tree["coeff"][i][4], 0.0)
    else:
        fit_y = legendre(fit_x, tree["coeff"][i][0], tree["coeff"][i][1], tree["coeff"][i][2], tree["coeff"][i][3], tree["coeff"][i][4], tree["coeff"][i][5])

    ax.plot(fit_x, fit_y, "C0", ls = "dashed")

    plt.subplots_adjust(left = 0.15, right=0.99, top=0.93, bottom = 0.13, hspace=0.01, wspace=0.01)
    ax.set_title(arg_dict["title"], fontsize = 40)
    ax.set_xlabel(r"$\cos \theta_{\rm CM}$", fontsize = 40)
    ax.set_ylabel(r"d$\sigma$/d$\Omega$ [mb/sr]", fontsize = 40)
    
    img_path = os.path.join(script_dir, "../results/img/yield/{}.{}".format(arg_dict["savename"], arg_dict["img_type"]))
    plt.savefig(img_path, format=arg_dict["img_type"], dpi=600, transparent=True)
    
    plt.show()


if __name__ == '__main__':

    img_type = "png"

    # +----+
    # | Kp |
    # +----+
    arg_dict = {
        "rootdata" : "for_cusp/dcs_kp_9999_for_cusp.root",
        "idx" : 54,
        "title": r"$K^-p\rightarrow K^-p$",
        "savename": "kp_legendre_result_example",
        "img_type": img_type,
        "coeff_flag": False
    }
    plot(arg_dict)

    # +-----+
    # | K0n |
    # +-----+
    arg_dict = {
        "rootdata" : "for_cusp/dcs_k0n_211_for_cusp.root",
        "idx" : 55,
        "title": r"$K^-p \rightarrow \overline{K}^0n$",
        "savename": "k0n_legendre_result_example",
        "img_type": img_type,
        "coeff_flag": False
    }
    plot(arg_dict)

    # +---------------+
    # | Kp->pi+Sigma- |
    # +---------------+
    arg_dict = {
        "rootdata" : "for_cusp/dcs_pi+sigma-_2112_for_cusp.root",
        "idx" : 55,
        "title": r"$K^-p \rightarrow \pi^+\Sigma^-$",
        "savename": "pi+sigma-_legendre_result_example",    
        "img_type": img_type,
        "coeff_flag": True
    }
    plot(arg_dict)