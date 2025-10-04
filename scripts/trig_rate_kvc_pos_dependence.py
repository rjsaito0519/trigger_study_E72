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

def f(x, amp, tau, c):
    return amp*np.exp( x/tau ) + c

dir = ""

data = []
momenta = np.asarray([645.0, 665.0, 685.0, 705.0, 725.0, 745.0, 765.0, 790.0, 814.0, 842.0, 870.0, 902.0, 933.0])

# for mom in list(range(605, 1006, 20)):
for mom in momenta:
    root_file_path = os.path.join(script_dir, f"../data/kvc_pos_optimize/{dir}/beam_mom{mom:.0f}.root")
    if os.path.exists(root_file_path):
        file = uproot.open(root_file_path)
        tree = file["tree"].arrays(library="np")
        print(tree["n_trig"])
        data.append([
            mom,
            tree["kaon_rate2024"][0]*tree["n_trig_all"][0]/tree["n_kaon_all"][0],
            tree["kaon_rate2024"][0]*tree["n_trig"][0][0]/tree["n_kaon_all"][0], # no decay
            tree["kaon_rate2024"][0]*tree["n_trig"][0][1]/tree["n_kaon_all"][0], # muon
            tree["kaon_rate2024"][0]*(tree["n_trig"][0][2]+tree["n_trig"][0][3])/tree["n_kaon_all"][0], # pion
        ])
data = np.asarray(data)


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
ax.plot(data[:, 0], data[:, 1], "o", color = "C0")

# # --------------------
# model = lfm.Model(f, name = "f")
# params = model.make_params()
# params.add('amp', value=2.74e4, min=0)
# # params.add('x0',  value=700)
# params.add('tau', value=115, min=1)
# params.add('c',  value=1.0)
# result = model.fit(x = data[:, 0], data = data[:, 1], params=params, method='leastsq')
# print(result.fit_report())
# fit_x = np.linspace(600, 1000, 10000)
# fit_y = result.eval_components(x=fit_x)["f"]
# ax.plot(fit_x, fit_y, "--")

# target = 1.5 # kHz
# max_iter = 50
# tol = 1e-12,
# n_step = 100
# min_mom = data[:, 0].min()
# max_mom = data[:, 0].max()
# target_mom = None
# for _ in range(max_iter):
#     moms = np.linspace(min_mom, max_mom, n_step)
#     rate = result.eval_components(x=moms)["f"]
#     diff = np.abs( np.full_like(rate, target) - rate )
#     target_mom = moms[np.argmin(diff)]
#     mom_range = max_mom - min_mom
#     min_mom = target_mom - mom_range/n_step * 2
#     max_mom = target_mom + mom_range/n_step * 2
#     if np.min(diff) < tol:
#         break
# print(target_mom, result.eval_components(x=[884.5])["f"])
# # --------------------


# plt.legend(fontsize = 18)
ax.set_xlabel(r"$p_K$ [MeV/c]")
ax.set_ylabel("Accidental Trigger Rate [kHz]")
# img_save_path = os.path.join(script_dir, "../results/img/01.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()



fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
ax.plot(data[:, 0], data[:, 2], "--o", color = "C0", label = "no decay")
ax.plot(data[:, 0], data[:, 3], "--s", color = "C1", label = r"decay to $\mu$")
ax.plot(data[:, 0], data[:, 4], "--^", color = "C2", label = r"decay to $\pi$")
# plt.legend(fontsize = 18)
ax.set_xlabel(r"$p_K$ [MeV/c]")
ax.set_ylabel("Accidental Trigger Rate [kHz]")
# img_save_path = os.path.join(script_dir, "../results/img/02.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


sys.exit()

data = {
    645: [],
    735: [],
    805: []
}

for mom in [645, 735, 805]:
    for z_pos in range(45, 101, 5):
        root_file_path = os.path.join(script_dir, f"../data/kvc_pos_optimize/{dir}/kvc_pos_optimize_beam_mom{mom}_{z_pos}.root")
        if os.path.exists(root_file_path):
            file = uproot.open(root_file_path)
            tree = file["tree"].arrays(library="np")
            print(tree["n_trig"])
            data[mom].append([
                z_pos,
                tree["kaon_rate2024"][0]*tree["n_trig_all"][0]/tree["n_kaon_all"][0],
                tree["kaon_rate2024"][0]*tree["n_trig"][0][0]/tree["n_kaon_all"][0], # no decay
                tree["kaon_rate2024"][0]*tree["n_trig"][0][1]/tree["n_kaon_all"][0], # muon
                tree["kaon_rate2024"][0]*(tree["n_trig"][0][2]+tree["n_trig"][0][3])/tree["n_kaon_all"][0], # pion
            ])
    data[mom] = np.array(data[mom])


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
for i, mom in enumerate([645, 735, 805]):
    ax.plot(data[mom][:, 0], data[mom][:, 1], "o", color = f"C{i}", label = r"$p_K = ${:.0f}".format(mom))
# plt.legend(fontsize = 18)
ax.set_xlabel("KVC Z position[mm]")
ax.set_ylabel("Accidental Trigger Rate [kHz]")
img_save_path = os.path.join(script_dir, "../results/img/01.png")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()

fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
for i, mom in enumerate([645, 735, 805]):
    if mom == 645:
        ax.plot(data[mom][:, 0], data[mom][:, 2], "--o", color = f"C{i}", label = "no decay")
        ax.plot(data[mom][:, 0], data[mom][:, 3], "--s", color = f"C{i}", label = r"decay to $\mu$")
        ax.plot(data[mom][:, 0], data[mom][:, 4], "--^", color = f"C{i}", label = r"decay to $\pi$")
    else:
        ax.plot(data[mom][:, 0], data[mom][:, 2], "--o", color = f"C{i}")
        ax.plot(data[mom][:, 0], data[mom][:, 3], "--s", color = f"C{i}")
        ax.plot(data[mom][:, 0], data[mom][:, 4], "--^", color = f"C{i}")
# plt.legend(fontsize = 18)
ax.set_xlabel("KVC Z position[mm]")
ax.set_ylabel("Accidental Trigger Rate [kHz]")
img_save_path = os.path.join(script_dir, "../results/img/02.png")
os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


el735 = np.array([
    [45, 740596, 700895, 94.6393, 67.3715, 93.6230],
    [50, 740230, 703503, 95.0384, 67.7456, 94.0117],
    [55, 741529, 707043, 95.3493, 67.8578, 94.3333],
    [60, 741781, 708899, 95.5672, 68.0852, 94.5612],
    [65, 742186, 710845, 95.7772, 68.3424, 94.7684],
    [70, 741819, 712131, 95.9979, 68.4957, 94.9598],
    [75, 742827, 714393, 96.1722, 68.5631, 95.1659],
    [80, 742049, 714350, 96.2672, 68.6681, 95.2455],
    [85, 742619, 716129, 96.4329, 68.7677, 95.4064],
    [90, 741951, 716361, 96.5510, 68.8683, 95.5248],
    [95, 740999, 716076, 96.6366, 68.9932, 95.5835],
    [100,742782, 718530, 96.7350, 68.9733, 95.7109]
])

# etaLambda 805
el805 = np.array([
    [45, 1115716, 1027754, 92.1161, 66.8629, 89.4558],
    [50, 1117171, 1034022, 92.5572, 67.3819, 89.8768],
    [55, 1114809, 1034781, 92.8214, 67.5409, 90.1443],
    [60, 1116626, 1040150, 93.1512, 67.8851, 90.4637],
    [65, 1114396, 1040208, 93.3428, 68.0972, 90.6343],
    [70, 1115347, 1043204, 93.5318, 68.2497, 90.8109],
    [75, 1115520, 1044468, 93.6306, 68.2202, 90.9016],
    [80, 1114409, 1044747, 93.7490, 68.3639, 91.0253],
    [85, 1116545, 1047308, 93.7990, 68.5197, 91.0768]
])

fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
ax.plot(el735[:, 0]*10, el735[:, 3], "--s", color = "C1")
ax.plot(el805[:, 0]*10, el805[:, 3], "--^", color = "C2")

# plt.legend(fontsize = 18)
ax.set_xlabel("KVC Z position[mm]")
ax.set_ylabel("Acceptance [%]")
plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.12)
# img_save_path = os.path.join(script_dir, "../results/img/02.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


# kp 645
kp645 = np.array([
# "Z", "n_entry", "n_trig", "acceptance", "acceptance_mp2", "acceptance_htofp"
    [45, 1616343, 1185908, 73.3698, 53.6799, 67.4502],
    [50, 1617095, 1186994, 73.4029, 53.7148, 67.4284],
    [55, 1616597, 1186525, 73.3965, 53.7137, 67.4250],
    [60, 1616113, 1185992, 73.3855, 53.6982, 67.4567],
    [65, 1617076, 1186306, 73.3612, 53.6726, 67.3391],
    [70, 1615997, 1186849, 73.4438, 53.7383, 67.4556],
    [75, 1616925, 1187138, 73.4195, 53.7821, 67.4355],
    [80, 1617569, 1187576, 73.4173, 53.7553, 67.4410],
    [85, 1617510, 1187540, 73.4178, 53.7687, 67.4196],
    [90, 1616148, 1187634, 73.4855, 53.7532, 67.5038],
    [95, 1617132, 1187677, 73.4434, 53.7611, 67.4450],
    [100,1617209, 1187843, 73.4502, 53.7430, 67.4204]
])

# kp 735
kp735 = np.array([
    [45, 1694483, 1221273, 72.0735, 55.4891, 64.2804],
    [50, 1695983, 1224544, 72.2026, 55.5471, 64.3511],
    [55, 1695736, 1225149, 72.2488, 55.5929, 64.3873],
    [60, 1696123, 1226461, 72.3097, 55.6493, 64.4053],
    [65, 1695008, 1225148, 72.2798, 55.6378, 64.4425],
    [70, 1694139, 1225342, 72.3283, 55.6762, 64.4240],
    [75, 1695452, 1225822, 72.3006, 55.6737, 64.4213],
    [80, 1695466, 1225949, 72.3075, 55.6625, 64.4398],
    [85, 1695364, 1226211, 72.3273, 55.6721, 64.4476],
    [90, 1694510, 1225257, 72.3075, 55.7161, 64.4061],
    [95, 1696120, 1227428, 72.3668, 55.6942, 64.4719],
    [100,1695985, 1226861, 72.3391, 55.7271, 64.4634]
])

# kp 805
kp805 = np.array([
    [45, 1741495, 1213453, 69.6788, 54.8935, 59.6030],
    [50, 1741251, 1218178, 69.9599, 55.1576, 59.7810],
    [55, 1740681, 1222304, 70.2199, 55.4394, 59.8355],
    [60, 1740480, 1226729, 70.4822, 55.6474, 59.9393],
    [65, 1740187, 1228469, 70.5941, 55.7612, 59.9913],
    [70, 1740546, 1231442, 70.7503, 55.9277, 60.0302],
    [75, 1740627, 1234046, 70.8966, 56.0356, 60.0953],
    [80, 1741826, 1236813, 71.0067, 56.1891, 60.1554]
])

fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
ax.plot(kp645[:, 0]*10, kp645[:, 3], "--o", color = "C0")
ax.plot(kp735[:, 0]*10, kp735[:, 3], "--s", color = "C1")
ax.plot(kp805[:, 0]*10, kp805[:, 3], "--^", color = "C2")

# plt.legend(fontsize = 18)
ax.set_xlabel("KVC Z position[mm]")
ax.set_ylabel("Acceptance [%]")
plt.subplots_adjust(left = 0.15, right=0.98, top=0.98, bottom = 0.12)
# img_save_path = os.path.join(script_dir, "../results/img/02.png")
# os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
# plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()