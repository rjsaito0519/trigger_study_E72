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

color_dict = dict(
    bubble       = "C0",
    bubble1970   = "C0",
    bubble1968   = "C0",
    bubble1975   = "C0",
    deuteron_bubble1 = "C4",
    deuteron_bubble2 = "C4",
    crystal_ball = "C1",
    counter      = "C2",
    spark        = "C6",
)
marker_dict = dict(
    bubble       = "o",
    bubble1970   = "o",
    bubble1968   = "*",
    bubble1975   = "d",
    deuteron_bubble1 = "P",
    deuteron_bubble2 = "X",
    crystal_ball = "s",
    counter      = "^",
    spark        = "p",
)

def sqrt_s(mom_kaon):
    return np.sqrt( (np.sqrt(mom_kaon**2 +  np.full_like(mom_kaon, mass_kaon)**2) + np.full_like(mom_kaon, mass_proton))**2 - mom_kaon**2 )

def sqrt_s_err(mom_kaon, mom_kaon_err):
    numerator = (2 * mom_kaon * (np.sqrt(mass_kaon**2 + mom_kaon**2) + mass_proton)) / np.sqrt(mass_kaon**2 + mom_kaon**2) - 2 * mom_kaon
    denominator = 2 * np.sqrt((np.sqrt(mass_kaon**2 + mom_kaon**2) + mass_proton)**2 - mom_kaon**2)
    derivative = numerator / denominator
    return np.sqrt( derivative**2 * mom_kaon_err**2 )


def plot(arg_dict1, arg_dict2, arg_dict3):

    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    
    plt.subplots_adjust(left = 0.09, right=0.99, top=0.9, bottom = 0.1)

    for ax, ax_sub, arg_dict in zip([ax1, ax2, ax3], [ax4, ax5, ax6], [arg_dict1, arg_dict2, arg_dict3]):
        root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["rootdata"]))
        file = uproot.open(root_file_path)
        tree = file["tree"].arrays(library="np")

        for rawdata, key in zip(arg_dict["data"], arg_dict["data_key"]):
            mask = (range_left < sqrt_s(rawdata[:, 0])) * (sqrt_s(rawdata[:, 0]) < range_right)
            data = rawdata[mask]
            ax.errorbar(
                sqrt_s(data[:, 0])/1000, data[:, 2], 
                xerr = sqrt_s_err(data[:, 0], data[:, 1])/1000, 
                yerr = data[:, 3], 
                fmt = marker_dict[key], 
                capsize = 0, 
                markeredgecolor = color_dict[key], 
                ms = 12 if marker_dict[key] == "*" else 8, 
                ecolor='gray', 
                color="w", 
                markeredgewidth = 1.5, 
                label = key,
                zorder = 3
            )
            
            # -- E72 -----
            a0     = [item[0] for item in tree["coeff"]]
            a0_err = [item[0] for item in tree["coeff_err"]]
            ax.errorbar(
                sqrt_s(tree["mom"])/1000, a0, 
                xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
                yerr = a0_err, 
                fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
            )
        
            if key == "bubble1970":
                # -- bubble -----
                ax_sub.errorbar(
                    sqrt_s(data[:, 0])/1000, data[:, 2], 
                    xerr = sqrt_s_err(data[:, 0], data[:, 1])/1000, 
                    yerr = data[:, 3], 
                    fmt = marker_dict[key], 
                    capsize = 0, 
                    markeredgecolor = color_dict[key], 
                    ms = 12 if marker_dict[key] == "*" else 8, 
                    ecolor='gray', 
                    color="w", 
                    markeredgewidth = 1.5, 
                    label = key,
                    zorder = 3
                )

                # -- E72 -----
                ax_sub.errorbar(
                    sqrt_s(tree["mom"])/1000, a0, 
                    xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
                    yerr = a0_err, 
                    fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
                )
            
        ax.axvline(mass_tot/1000, ls = "dashed", color="red", lw = 1.0, zorder = 0)
        ax.grid(color="gray", linestyle="dotted", linewidth=1)
        offset = 0.15*0.1
        ax.set_xlim(range_left/1000-offset, range_right/1000+offset)

        ax_sub.axvline(mass_tot/1000, ls = "dashed", color="red", lw = 1.0, zorder = 0)
        ax_sub.grid(color="gray", linestyle="dotted", linewidth=1)
        xmin = 1.63607743
        xmax = 1.68912484
        offset = (xmax-xmin)*0.15
        ax_sub.set_xlim(xmin-offset, xmax+offset)

        up_left = fig.transFigure.inverted().transform(
            ax.transData.transform(
                (xmin-offset, ax.get_ylim()[0])
            )
        )
        up_right = fig.transFigure.inverted().transform(
            ax.transData.transform(
                (xmax+offset, ax.get_ylim()[0])
            )
        )
        
        low_left = fig.transFigure.inverted().transform(
            ax_sub.transData.transform(
                (xmin-offset, ax_sub.get_ylim()[1])
            )
        )
        low_right = fig.transFigure.inverted().transform(
            ax_sub.transData.transform(
                (xmax+offset, ax_sub.get_ylim()[1])
            )
        )

        line_props = dict(color="k", linestyle="--", linewidth=1.5, zorder=0)
        fig.add_artist(plt.Line2D([up_left[0], low_left[0]], [up_left[1], low_left[1]], **line_props))  # 左の線
        fig.add_artist(plt.Line2D([up_right[0], low_right[0]], [up_right[1], low_right[1]], **line_props))  # 左の線
        

    ax1.set_title(r"$K^-p\rightarrow K^-p$")
    ax2.set_title(r"$K^-p\rightarrow \overline{K}^0n$")
    ax3.set_title(r"$K^-p\rightarrow \pi^+\Sigma^-$")
    fig.supxlabel(r"$W$ [GeV]", x = 0.54)
    ax1.set_ylabel(r"$A_0$")
    ax4.set_ylabel(r"$A_0$")

    img_path = os.path.join(script_dir, "../results/img/yield/Expected_A0_plot.pdf")
    plt.savefig(img_path, format='pdf', dpi=600, transparent=True)
    img_path = os.path.join(script_dir, "../results/img/yield/Expected_A0_plot.png")
    plt.savefig(img_path, format='png', dpi=600, transparent=True)

    plt.show()

def plot_up_panel(arg_dict1, arg_dict2, arg_dict3):

    fig = plt.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    plt.subplots_adjust(left = 0.09, right=0.99, top=0.88, bottom = 0.17)

    for ax, arg_dict in zip([ax1, ax2, ax3], [arg_dict1, arg_dict2, arg_dict3]):
        root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["rootdata"]))
        file = uproot.open(root_file_path)
        tree = file["tree"].arrays(library="np")

        for rawdata, key in zip(arg_dict["data"], arg_dict["data_key"]):
            mask = (range_left < sqrt_s(rawdata[:, 0])) * (sqrt_s(rawdata[:, 0]) < range_right)
            data = rawdata[mask]
            ax.errorbar(
                sqrt_s(data[:, 0])/1000, data[:, 2], 
                xerr = sqrt_s_err(data[:, 0], data[:, 1])/1000, 
                yerr = data[:, 3], 
                fmt = marker_dict[key], 
                capsize = 0, 
                markeredgecolor = color_dict[key], 
                ms = 12 if marker_dict[key] == "*" else 8, 
                ecolor='gray', 
                color="w", 
                markeredgewidth = 1.5, 
                label = key,
                zorder = 3
            )
            
            # -- E72 -----
            a0     = [item[0] for item in tree["coeff"]]
            a0_err = [item[0] for item in tree["coeff_err"]]
            ax.errorbar(
                sqrt_s(tree["mom"])/1000, a0, 
                xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
                yerr = a0_err, 
                fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
            )
                    
        ax.axvline(mass_tot/1000, ls = "dashed", color="red", lw = 1.0, zorder = 0)
        ax.grid(color="gray", linestyle="dotted", linewidth=1)
        offset = 0.15*0.1
        ax.set_xlim(range_left/1000-offset, range_right/1000+offset)
        
    ax1.set_title(r"$K^-p\rightarrow K^-p$")
    ax2.set_title(r"$K^-p\rightarrow \overline{K}^0n$")
    ax3.set_title(r"$K^-p\rightarrow \pi^+\Sigma^-$")
    fig.supxlabel(r"$W$ [GeV]", x = 0.54)
    ax1.set_ylabel(r"$A_0$", fontsize = 35)

    img_path = os.path.join(script_dir, "../results/img/yield/Expected_A0_plot_uppanel.pdf")
    plt.savefig(img_path, format='pdf', dpi=600, transparent=True)
    img_path = os.path.join(script_dir, "../results/img/yield/Expected_A0_plot_uppanel.png")
    plt.savefig(img_path, format='png', dpi=600, transparent=True)
    plt.show()


def legend_plot():
    fig = plt.figure(figsize=(5, 5))
    ax   = fig.add_subplot(111)
    kp_spark         = np.genfromtxt(os.path.join(script_dir,"../data/legendre/Kp_spark_chamber_legendre.csv"), skip_header=1, delimiter=",")
    k0n_bubble1      = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_bubble_chamber1970_legendre.csv"), skip_header=1, delimiter=",")
    k0n_bubble2      = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_bubble_chamber1968_legendre.csv"), skip_header=1, delimiter=",")
    k0n_bubble3      = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_bubble_chamber1975_legendre.csv"), skip_header=1, delimiter=",")
    k0n_crystal_ball = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_crystal_ball_legendre.csv"), skip_header=1, delimiter=",")
    k0n_counter      = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_counter_legendre.csv"), skip_header=1, delimiter=",")

    root_file_path = os.path.join(script_dir, "../results/root/dcs_kp_9999.root")
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    legends = []
    i = 0
    for data, key in zip([k0n_bubble1, k0n_bubble2, k0n_bubble3, k0n_crystal_ball, k0n_counter, kp_spark], ["bubble1970", "bubble1968", "bubble1975", "crystal_ball", "counter", "spark"]):
        buf = ax.errorbar(sqrt_s(data[:, 0]), data[:, 2*i+2], xerr = sqrt_s_err(data[:, 0], data[:, 1]), yerr = data[:, 2*i+3], fmt = marker_dict[key], capsize = 0, markeredgecolor = color_dict[key], ms = 12, ecolor='k', color="w", markeredgewidth = 1.5, label = key)
        legends.append(buf)

    # -- E72 -----
    a0     = [item[0] for item in tree["coeff"]]
    a0_err = [item[0] for item in tree["coeff_err"]]
    buf = ax.errorbar(
        sqrt_s(tree["mom"])/1000, a0, 
        xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
        yerr = a0_err, 
        fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
    )
    legends.append(buf)

    fig_legend = plt.figure(figsize=(8, 8))  # 凡例専用のFigureを作成
    fig_legend.legend(
        handles=legends,
        labels=["bubble chamber (1970)", "bubble chamber (1968)", "bubble chamber (1975)", "Crystal Ball", "counter", "spark chamber", "J-PARC E72 (expected)"],   # ラベル（凡例テキスト）
        loc='center', handletextpad = 0.5, handlelength=0.5
    )    
    fig_legend.savefig(os.path.join(script_dir, "../results/img/yield/A0_legend.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
    fig_legend.savefig(os.path.join(script_dir, "../results/img/yield/A0_legend.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)

    plt.close(fig_legend)
    plt.close(fig)

def plot_only_kp(arg_dict):
    fig = plt.figure(figsize=(5, 5))
    ax   = fig.add_subplot(111)
    kp_spark         = np.genfromtxt(os.path.join(script_dir,"../data/legendre/Kp_spark_chamber_legendre.csv"), skip_header=1, delimiter=",")
    k0n_bubble1      = np.genfromtxt(os.path.join(script_dir,"../data/legendre/K0n_bubble_chamber1970_legendre.csv"), skip_header=1, delimiter=",")

    root_file_path = os.path.join(script_dir, "../results/root/dcs_kp_9999.root")
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    legends = []
    i = 0
    for data, key in zip([k0n_bubble1, kp_spark], ["bubble1970", "spark"]):
        buf = ax.errorbar(sqrt_s(data[:, 0]), data[:, 2*i+2], xerr = sqrt_s_err(data[:, 0], data[:, 1]), yerr = data[:, 2*i+3], fmt = marker_dict[key], capsize = 0, markeredgecolor = color_dict[key], ms = 12, ecolor='k', color="w", markeredgewidth = 1.5, label = key)
        legends.append(buf)

    # -- E72 -----
    a0     = [item[0] for item in tree["coeff"]]
    a0_err = [item[0] for item in tree["coeff_err"]]
    buf = ax.errorbar(
        sqrt_s(tree["mom"])/1000, a0, 
        xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
        yerr = a0_err, 
        fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
    )
    legends.append(buf)

    fig_legend = plt.figure(figsize=(8, 8))  # 凡例専用のFigureを作成
    fig_legend.legend(
        handles=legends,
        labels=["bubble chamber", "spark chamber", "J-PARC E72 (expected)"],   # ラベル（凡例テキスト）
        loc='center', handletextpad = 0.5, handlelength=0.5
    )    
    fig_legend.savefig(os.path.join(script_dir, "../results/img/yield/kp_A0_legend.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)

    plt.close(fig_legend)
    plt.close(fig)


    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111)

    root_file_path = os.path.join(script_dir, "../results/root/{}".format(arg_dict["rootdata"]))
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    for rawdata, key in zip(arg_dict["data"], arg_dict["data_key"]):
        mask = (range_left < sqrt_s(rawdata[:, 0])) * (sqrt_s(rawdata[:, 0]) < range_right)
        data = rawdata[mask]
        ax.errorbar(
            sqrt_s(data[:, 0])/1000, data[:, 2], 
            xerr = sqrt_s_err(data[:, 0], data[:, 1])/1000, 
            yerr = data[:, 3], 
            fmt = marker_dict[key], 
            capsize = 0, 
            markeredgecolor = color_dict[key], 
            ms = 12 if marker_dict[key] == "*" else 8, 
            ecolor='gray', 
            color="w", 
            markeredgewidth = 1.5, 
            label = key,
            zorder = 3
        )
        
        # -- E72 -----
        a0     = [item[0] for item in tree["coeff"]]
        a0_err = [item[0] for item in tree["coeff_err"]]
        ax.errorbar(
            sqrt_s(tree["mom"])/1000, a0, 
            xerr = sqrt_s_err(tree["mom"], tree["mom_err"])/1000, 
            yerr = a0_err, 
            fmt = "o", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color = "w", zorder = 5, elinewidth = 2, label = "E72"
        )
                
    ax.axvline(mass_tot/1000, ls = "dashed", color="red", lw = 1.5, zorder = 0)
    ax.grid(color="gray", linestyle="dotted", linewidth=1)
    offset = 0.15*0.1
    # ax.set_xlim(range_left/1000-offset, range_right/1000+offset)
    ax.set_xlim(mass_tot/1000 - 0.075, mass_tot/1000 + 0.075)
        
    plt.subplots_adjust(left = 0.15, right=0.95, top=0.88, bottom = 0.17)
    ax.set_title(r"$K^-p\rightarrow K^-p$")
    ax.set_xlabel(r"$W$ [GeV]")
    ax.set_ylabel(r"$A_0$", fontsize = 35)

    img_path = os.path.join(script_dir, "../results/img/yield/kp_Expected_A0_plot_uppanel.png")
    plt.savefig(img_path, format='png', dpi=600, transparent=True)
    plt.show()


if __name__ == '__main__':

    legend_plot()

    # +----+
    # | Kp |
    # +----+
    arg_dict1 = {
        "data": [
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/Kp_bubble_chamber1970_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/Kp_spark_chamber_legendre.csv"), skip_header=1, delimiter=","),
        ],
        "rootdata": "dcs_kp_9999_for_cusp_for_cusp.root",
        "data_key": ["bubble1970", "spark"],
        "title": r"$K^-p\rightarrow K^-p$",
    }

    # plot_only_kp(arg_dict1)
    # sys.exit()

    # +-----+
    # | K0n |
    # +-----+
    arg_dict2 = {
        "data": [
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/K0n_bubble_chamber1970_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/K0n_bubble_chamber1968_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/K0n_bubble_chamber1975_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/K0n_crystal_ball_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/K0n_counter_legendre.csv"), skip_header=1, delimiter=","),
            
        ],
        "rootdata": "dcs_k0n_211_for_cusp_for_cusp.root",
        "data_key": ["bubble1970", "bubble1968", "bubble1975", "crystal_ball", "counter"],
        "title": r"$K^-p \rightarrow \overline{K}^0n$",
        "savename": "k0n_coeff_result"
    }

    # +---------------+
    # | Kp->pi+Sigma- |
    # +---------------+
    arg_dict3 = {
        "data": [
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/pi+Sigma-_bubble_chamber1970_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/pi+Sigma-_bubble_chamber1968_legendre.csv"), skip_header=1, delimiter=","),
            np.genfromtxt(os.path.join(script_dir, "../data/legendre/pi+Sigma-_bubble_chamber1975_legendre.csv"), skip_header=1, delimiter=","),            
        ],
        "rootdata": "dcs_pi+sigma-_2112_for_cusp_for_cusp.root",
        "data_key": ["bubble1970", "bubble1968", "bubble1975"],
        "title": r"$K^-p \rightarrow \pi^+\Sigma^-$",
        "savename": "pi+sigma-_coeff_result"
    }
    # plot(arg_dict1, arg_dict2, arg_dict3)
    plot_up_panel(arg_dict1, arg_dict2, arg_dict3)
    