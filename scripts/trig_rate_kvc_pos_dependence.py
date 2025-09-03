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

mom = 735

for z_pos in range(45, 81, 5):
    root_file_path = os.path.join(script_dir, f"../data/kvc_pos_optimize/kvc_pos_optimize_beam_mom{mom}_{z_pos}.root")
    if os.path.exists(root_file_path):
        file = uproot.open(root_file_path)
        tree = file["tree"].arrays(library="np")

        print(tree["n_trig_all"]/tree["n_kaon_all"])