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

data = np.array([
    [700, 27.4*10**3],
    [735, 38.4*10**3],
    [750, 44.0*10**3]
])

model = lfm.LinearModel()
params = model.guess(x = data[:, 0], data = data[:, 1] )
result = model.fit(x = data[:, 0], data = data[:, 1], params=params, method='leastsq')
print(result.fit_report())
fit_x = np.linspace(650, 800, 2)
fit_y = result.eval_components(x=fit_x)["linear"]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(fit_x, fit_y, color = "C1")
ax.plot(data[:, 0], data[:, 1], "o", ms = 10, color = "C0")

ax.set_xlabel(r"$p_K$ [MeV/c]")
ax.set_ylabel(r"$K^-$ Beam Intensity [$\times 10^3$/spill]")
ax.yaxis.set_major_formatter(ptick.EngFormatter())

plt.subplots_adjust(left = 0.18, right=0.98, top=0.98, bottom = 0.12)

img_save_path = os.path.join( script_dir, "../results/img/yield/beam_intensity.pdf")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)

plt.show()

