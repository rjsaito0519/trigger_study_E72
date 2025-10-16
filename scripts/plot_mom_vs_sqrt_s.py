import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import numpy as np
import uproot
import os
import sys

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

mass_Lambda = 1115.683
mass_eta    =  547.862
mass_kaon   =  493.677
mass_proton =  938.27208816

mass_tot = mass_Lambda + mass_eta

def f_sqrt_s(mom_kaon):
    return np.sqrt( (np.sqrt(mom_kaon**2 +  np.full_like(mom_kaon, mass_kaon)**2) + np.full_like(mom_kaon, mass_proton))**2 - mom_kaon**2 )


fig = plt.figure(figsize=(8, 6))
ax  = fig.add_subplot(111)

mom = np.linspace(600, 1000, 10000)
sqrt_s = f_sqrt_s(mom)

ax.plot(mom, sqrt_s)
ax.axvline(620, ls = "dashed", color = "red")
ax.axvline(955, ls = "dashed", color = "red")

momenta = np.asarray([645.0, 665.0, 685.0, 705.0, 725.0, 745.0, 765.0, 790.0, 814.0, 842.0, 870.0, 902.0, 933.0, 735.0])
ax.plot(momenta, f_sqrt_s(momenta), "o", color = "C1")

ax.set_xlabel(r"$p_K$ [MeV/c]")
ax.set_ylabel(r"$\sqrt{s}$ [MeV]")

plt.subplots_adjust(left = 0.18, right=0.98, top=0.98, bottom = 0.16)
img_save_path = os.path.join(script_dir, "../results/img/mom_vs_sqrt_s.png")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)

plt.show()