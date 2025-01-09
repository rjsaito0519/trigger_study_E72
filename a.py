import matplotlib.pyplot as plt
import numpy as np

# 7x7のsubplotを作成
fig, axes = plt.subplots(10, 6, figsize=(10, 10))

# 各subplotにデータをプロット
for i, ax in enumerate(axes.flat):  # axes.flatで1次元のリストとして扱う
    # データを生成 (例としてランダムなsin波)
    x = np.linspace(0, 10, 100)
    y = np.sin(x + i)
    ax.plot(x, y)
    ax.set_title(f"Plot {i+1}")
    ax.set_xticks([])  # X軸ラベルを削除（必要に応じて調整）
    ax.set_yticks([])  # Y軸ラベルを削除（必要に応じて調整）

# レイアウト調整
plt.tight_layout()
plt.show()
