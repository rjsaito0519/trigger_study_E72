#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time allocation optimizer:
- Equalize counts r_i * t_i under total time T
- Bounds: lower <= t_i <= upper
- Quantize to fixed step (e.g., 10 minutes) while preserving sum/bounds

Example matches your case:
  T = 4.5 days * 24 h = 108 h
  momenta = 645, 665, ..., 925 (step 20)
  rate(x) = 137.962111*exp(x/126.840771) - 6979.77790   [/spill]
  bounds: lower=1 h, upper=12 h
  step: 10 min (1/6 h)
"""

from __future__ import annotations
import math
import numpy as np
from dataclasses import dataclass
from typing import Tuple, List, Optional

# ---------- Core water-filling (continuous) ----------
def water_filling_with_bounds(
    r: np.ndarray,
    T: float,
    lower: float,
    upper: float,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    Solve for t_i = clip(K / r_i, lower, upper) so that sum t_i = T.
    Returns (t_cont, K_star, status), where status in {"at_lower","active","at_upper"}.
    """
    r = np.asarray(r, dtype=float)
    m = len(r)
    if np.any(r <= 0):
        raise ValueError("All rates r_i must be > 0.")

    min_sum, max_sum = m * lower, m * upper
    if not (min_sum - tol <= T <= max_sum + tol):
        raise ValueError(
            f"Total time T={T}h not achievable with bounds [{lower},{upper}] over {m} conditions: "
            f"min={min_sum}, max={max_sum}"
        )

    def total_time_for_K(K: float) -> Tuple[float, np.ndarray]:
        t = K / r
        t = np.clip(t, lower, upper)
        return float(t.sum()), t

    # Bracket K: ensure sum(K_hi) >= T
    K_lo = 0.0
    sum_lo, _ = total_time_for_K(K_lo)
    K_hi = max(1.0, upper * float(r.max()))
    sum_hi, _ = total_time_for_K(K_hi)
    it_grow = 0
    while sum_hi < T and it_grow < 60:
        K_hi *= 2.0
        sum_hi, _ = total_time_for_K(K_hi)
        it_grow += 1

    # Bisection
    for _ in range(max_iter):
        K_mid = 0.5 * (K_lo + K_hi)
        s_mid, _ = total_time_for_K(K_mid)
        if abs(s_mid - T) <= tol:
            K_star = K_mid
            break
        if s_mid > T:
            K_hi = K_mid
        else:
            K_lo = K_mid
    else:
        K_star = 0.5 * (K_lo + K_hi)

    _, t_cont = total_time_for_K(K_star)

    # Status flags
    status = np.where(
        t_cont <= lower + 1e-12, "at_lower",
        np.where(t_cont >= upper - 1e-12, "at_upper", "active")
    )
    return t_cont, K_star, status

# ---------- Quantization to fixed step preserving sum/bounds ----------
def quantize_to_step_preserve_sum(
    t_cont: np.ndarray,
    T: float,
    step: float,
    lower: float,
    upper: float,
) -> np.ndarray:
    """
    Largest-remainder method in units of `step`, respecting [lower, upper] and preserving sum to machine precision.
    """
    # Floor to step
    units_cont = t_cont / step
    units_floor = np.floor(units_cont)
    t_floor = units_floor * step
    t_floor = np.clip(t_floor, lower, upper)

    current_sum = float(t_floor.sum())
    diff_units = int(round((T - current_sum) / step))

    # Fractional parts for priority
    frac = units_cont - np.floor(units_cont)

    tq = t_floor.copy()

    # Add units: pick largest fractional parts among those below upper
    if diff_units > 0:
        candidates = np.where(tq < upper - 1e-12)[0]
        order = candidates[np.argsort(frac[candidates])[::-1]]  # largest frac first
        i = 0
        while diff_units > 0:
            if i >= len(order):
                # refresh candidates in case many hit the upper bound
                candidates = np.where(tq < upper - 1e-12)[0]
                if len(candidates) == 0:
                    break
                order = candidates[np.argsort(frac[candidates])[::-1]]
                i = 0
            k = order[i]
            if tq[k] + step <= upper + 1e-12:
                tq[k] += step
                diff_units -= 1
            i += 1

    # Remove units: pick smallest fractional parts among those above lower
    elif diff_units < 0:
        candidates = np.where(tq > lower + 1e-12)[0]
        order = candidates[np.argsort(frac[candidates])]  # smallest frac first
        i = 0
        while diff_units < 0:
            if i >= len(order):
                candidates = np.where(tq > lower + 1e-12)[0]
                if len(candidates) == 0:
                    break
                order = candidates[np.argsort(frac[candidates])]
                i = 0
            k = order[i]
            if tq[k] - step >= lower - 1e-12:
                tq[k] -= step
                diff_units += 1
            i += 1

    # Final tiny cleanup
    tq = np.clip(tq, lower, upper)
    # snap exactly to the step grid
    tq = np.round(tq / step) * step
    # Nudge if any tiny residual remains due to float
    residual_units = int(round((T - float(tq.sum())) / step))
    if residual_units != 0:
        if residual_units > 0:
            candidates = np.where(tq < upper - 1e-12)[0]
            order = candidates[np.argsort(frac[candidates])[::-1]]
            for k in order[:residual_units]:
                tq[k] += step
        else:
            candidates = np.where(tq > lower + 1e-12)[0]
            order = candidates[np.argsort(frac[candidates])]
            for k in order[:(-residual_units)]:
                tq[k] -= step
        tq = np.round(tq / step) * step
    return tq

# ---------- Helpers ----------
def hhmm(hours: float) -> str:
    h = int(math.floor(hours + 1e-9))
    m = int(round((hours - h) * 60))
    if m == 60:
        h += 1
        m = 0
    return f"{h:02d}:{m:02d}"

# ---------- Your concrete case ----------
def kaon_intensity(x: np.ndarray) -> np.ndarray:
    # [/spill], spills/hour factor cancels for inverse-rate allocation
    return np.clip(
        137.962111 * np.exp(x / 126.840771) - 6979.77790,
        None,
        137.962111 * np.exp(884.5 / 126.840771) - 6979.77790
        )

def main():
    # Problem setup
    days = 2.5
    T = days * 24.0  # 108 h
    lower, upper = 1.0, 12.0
    step = 10.0 / 60.0  # 10 minutes

    # momenta = np.arange(645, 926, 20, dtype=float)  # 645, 665, ..., 925
    momenta = np.asarray([645.0, 665.0, 685.0, 705.0, 725.0, 745.0, 765.0, 790.0, 814.0, 842.0, 870.0, 902.0, 933.0])
    r = kaon_intensity(momenta)
    if np.any(r <= 0):
        raise RuntimeError("Non-positive rate encountered; check function/range.")

    # Continuous solution (water-filling)
    t_cont, K_star, status = water_filling_with_bounds(r, T, lower, upper)

    # Quantize to 10-min units
    t_10 = quantize_to_step_preserve_sum(t_cont, T, step, lower, upper)

    # Sanity checks
    assert abs(t_10.sum() - T) < 1e-6
    assert np.all(t_10 >= lower - 1e-12)
    assert np.all(t_10 <= upper + 1e-12)
    assert np.allclose(t_10 / step, np.round(t_10 / step))

    # Report
    counts_cont = r * t_cont
    counts_10   = r * t_10
    avg_cont = counts_cont.mean()
    avg_10 = counts_10.mean()

    print(f"Total hours (target): {T:.6f} h")
    print(f"Sum continuous:      {t_cont.sum():.6f} h")
    print(f"Sum (10-min grid):   {t_10.sum():.6f} h")
    print(f"K* (continuous common count units): {K_star:.6f}")
    print()

    print(f"{'p [MeV/c]':>8}  {'rate(/spill)':>14}  {'t_cont[h]':>10}  {'t_10[h]':>8}  {'t_10[HH:MM]':>10}  {'status':>9}  {'counts_10(units)':>18}  {'Δ% from mean':>12}")
    for p, ri, tc, tq, st, c10 in zip(momenta, r, t_cont, t_10, status, counts_10):
        dperc = 100.0 * (c10 / avg_10 - 1.0)
        print(f"{int(p):8d}  {ri:14.3f}  {tc:10.3f}  {tq:8.3f}  {hhmm(tq):>10}  {st:>9}  {c10:18.3f}  {dperc:12.3f}")

    # Optional: write CSV
    # import pandas as pd
    # import pathlib
    # df = pd.DataFrame({
    #     "momentum_MeV_c": momenta.astype(int),
    #     "rate_per_spill": r,
    #     "t_cont_hours": t_cont,
    #     "t_10min_hours": t_10,
    #     "t_10min_hhmm": [hhmm(v) for v in t_10],
    #     "status": status,
    #     "counts_cont_units": counts_cont,
    #     "counts_10min_units": counts_10,
    # })
    # out = pathlib.Path("time_allocation_10min_bounds.csv")
    # df.to_csv(out, index=False)
    # print(f"\nSaved: {out.resolve()}")

if __name__ == "__main__":
    main()
