import matplotlib.pyplot as plt
import numpy as np
import uproot
import os
import glob
import sys
import statistics

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

hbarc = 197.3269804
mass_Lambda = 1115.683
mass_eta    =  547.862
mK = 493.677
mP = 938.27208816

cos_theta = np.array([
    [-0.95, 0.05],
    [-0.80, 0.10],
    [-0.60, 0.10],
    [-0.35, 0.15],
    [ 0.00, 0.20],
    [ 0.35, 0.15],
    [ 0.60, 0.10],
    [ 0.80, 0.10],
    [ 0.95, 0.05]
])

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

def legendre(x, a0, a1, a2, a3):
    term0 = 1
    term1 = x
    term2 = 3/2*x**2 - 1/2
    term3 = 5/2 * x**3 - 3/2 * x
    return a0*term0 + a1*term1 + a2*term2 + a3*term3



data = np.genfromtxt(os.path.join(script_dir, "../data/legendre/etaLambda_crystal_ball_legendre.csv"), skip_header=1, delimiter=",")
path_list = sorted(glob.glob(os.path.join(script_dir, "../data/dcs_data/kp_etalambda_mom*.csv")))

root_file_path = os.path.join(script_dir, "../results/root/dcs_eta_lambda_2212.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")


fig  = plt.figure(figsize=(10, 10))
ax01 = fig.add_subplot(5,3,1)
ax02 = fig.add_subplot(5,3,4)
ax03 = fig.add_subplot(5,3,7)
ax04 = fig.add_subplot(5,3,10)
ax05 = fig.add_subplot(5,3,13)
ax06 = fig.add_subplot(5,3,2)
ax07 = fig.add_subplot(5,3,5)
ax08 = fig.add_subplot(5,3,8)
ax09 = fig.add_subplot(5,3,11)
ax10 = fig.add_subplot(5,3,14)
ax11 = fig.add_subplot(5,3,3)
ax12 = fig.add_subplot(5,3,6)
ax13 = fig.add_subplot(5,3,9)
ax14 = fig.add_subplot(5,3,12)
ax15 = fig.add_subplot(5,3,15)
axs = [ax01, ax02, ax03, ax04, ax05, ax06, ax07, ax08, ax09, ax10, ax11, ax12, ax13, ax14]

i = 0
for ax, path in zip( axs, path_list ):
    # -- plot CB -----
    tmp_data = np.genfromtxt(path, skip_header=2, delimiter=",")
    mom = int(path[path.find("mom")+3:path.find(".csv")])

    CB = ax.errorbar(tmp_data[:, 0], tmp_data[:, 1], xerr = cos_theta[:, 1], yerr = tmp_data[:, 2], fmt = "s", capsize = 0, markeredgecolor = "red", ms = 0, ecolor='red', color='gray', markeredgewidth = 0.2, zorder = 3, label = "Crystal Ball")
    
    s = np.sqrt( (np.sqrt(mom**2 + mK**2) + mP)**2 - mom**2 )
    ax.text(0, 0.25, r"$p_K^{\rm lab}$" + " = {} MeV/c".format(mom), fontsize =18,  va = "center", ha = "center")
    ax.set_ylim(0, 0.31)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax)
    if ax not in [ax01, ax02, ax03, ax04, ax05]:
        ax.axes.yaxis.set_ticklabels([])
    if ax not in [ax05, ax10, ax14]:
        ax.axes.xaxis.set_ticklabels([])

    # -- plot E72 -----
    if mom != 770:
        index = np.where(tree["mom"] == mom)[0][0]
    else:
        index = np.where(tree["mom"] == 768)[0][0]

    E72 = ax.errorbar(
        tree["cos_theta"][index], tree["diff_cs"][index], xerr = tree["cos_theta_err"][index], yerr=tree["diff_cs_err"][index], 
        fmt = "s", capsize = 0, markeredgecolor = "k", ms = 0, ecolor='k', color = "w", zorder = 5, elinewidth = 1.5, label = "E72"
    )

    # -- for check error -----
    if mom == 734:
        print(tmp_data[:, 2])
        print(statistics.mean(tmp_data[:, 2]))
        print(tree["diff_cs_err"][index])
        print(statistics.mean(tree["diff_cs_err"][index]))

# -- for legend -----
ax15.legend(
    handles=[CB, E72],
    loc='lower center', handletextpad = 0.5, handlelength=0.5, fontsize = 24
)    
ax15.axis("off")

fig.suptitle(r"Differential Cross Section of $K^-p \rightarrow \eta\Lambda$", fontsize=26)
fig.supxlabel(r"$\cos \theta$")
fig.supylabel(r"d$\sigma$/d$\Omega$ [mb/sr]")


plt.subplots_adjust(left = 0.12, right=0.98, top=0.94, bottom = 0.1, hspace=0.02, wspace=0)
img_save_path = os.path.join(script_dir, "../results/img/yield/dcs_comparison_CB_and_E72.pdf")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

plt.show()