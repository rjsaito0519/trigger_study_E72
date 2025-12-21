import numpy as np
import math

def kaon_intensity(x: np.ndarray) -> np.ndarray:
    return np.clip(
        137.962111 * np.exp(x / 126.840771) - 6979.77790,
        1e-12,
        137.962111 * np.exp(884.5 / 126.840771) - 6979.77790
    )

def required_time_from_momentum(x, counts):
    """
    x      : momentum [MeV/c] (scalar or array-like)
    counts : desired event count per point

    return : required time [hours]
    """
    x = np.asarray(x, dtype=float)
    if counts <= 0:
        raise ValueError("counts must be positive")

    rate = kaon_intensity(x) / 4.24
    if np.any(rate <= 0):
        raise RuntimeError("non-positive rate encountered")

    return counts / rate


# ---------- 表示専用（1分丸め） ----------
def hhmm(hours: float) -> str:
    h = int(math.floor(hours + 1e-9))
    m = int(round((hours - h) * 60))
    if m == 60:
        h += 1
        m = 0
    return f"{h:02d}:{m:02d}"

data = [
    [645,	2.85E+07],
    [665,	5.48E+07],
    [685,	1.20E+08],
    [715,	4.05E+08],
    [735,	7.86E+08],
    [755,	3.67E+08],
    [790,	1.52E+08],
    [814,	1.57E+08],
    [842,	6.20E+07],
    [870,	4.70E+07],
    [933,   1.49E+08],
]

for it in data:
    t = required_time_from_momentum(it[0], it[1]) / 3600.0
    print(it[0], t)
