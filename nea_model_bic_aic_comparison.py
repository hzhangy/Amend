#!/usr/bin/env python3
"""
严格模型选择 v3：NEA平滑相变 vs MOND vs NFW（修正NFW公式）
NFW 使用标准公式：s = c*x，分子为 ln(1+s) - s/(1+s)
"""

import numpy as np
import glob
import os
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------
# 物理常数
# ----------------------------------------------------------------------
pc = 3.08567758e16       # m
km = 1000.0              # m
a0_mond = 1.2e-10        # m/s^2
H0_km_s_Mpc = 68.04
H0_km_s_kpc = H0_km_s_Mpc / 1000.0   # km/s / kpc

# ----------------------------------------------------------------------
# 数据加载
# ----------------------------------------------------------------------
files = sorted(glob.glob(os.path.join("Rotmod_LTG", "**", "*.dat"), recursive=True))
galaxies = []

for f in files:
    try:
        d = np.loadtxt(f, comments="#")
    except Exception:
        continue
    if d.ndim != 2 or d.shape[1] < 6:
        continue

    R = d[:, 0]
    Vobs = np.abs(d[:, 1])
    dVobs = np.abs(d[:, 2])
    Vgas = np.abs(d[:, 3])
    Vdisk = np.abs(d[:, 4])
    Vbulge = np.abs(d[:, 5]) if d.shape[1] >= 6 else 0.0

    mask = (R > 0) & (Vobs > 0) & (dVobs > 0) & (Vgas**2 + Vdisk**2 + Vbulge**2 > 0)
    if np.sum(mask) < 8:
        continue

    galaxies.append({
        "name": os.path.basename(f),
        "R": R[mask],
        "Vobs": Vobs[mask],
        "dVobs": dVobs[mask],
        "Vgas": Vgas[mask],
        "Vdisk": Vdisk[mask],
        "Vbulge": Vbulge[mask],
    })

print(f"Loaded usable galaxies: {len(galaxies)}")

# ----------------------------------------------------------------------
# 误差处理
# ----------------------------------------------------------------------
def get_eff_error(gal):
    return np.maximum(gal["dVobs"], 0.05 * gal["Vobs"])

# ----------------------------------------------------------------------
# 模型拟合函数
# ----------------------------------------------------------------------
def fit_nea_smooth(gal):
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon, Rc, n):
        Vbar2 = Vgas**2 + Upsilon * (Vdisk**2 + Vbulge**2)
        q = 1.0 + 1.0 / (1.0 + (R / Rc)**n)
        boost = (1.0 + R / Rc) ** (2.0 - q)
        return np.sqrt(np.maximum(Vbar2 * boost, 0.0))

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=[0.5, 5.0, 2.0],
                            bounds=([0.1, 0.5, 0.5], [1.5, 30.0, 6.0]),
                            maxfev=30000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        k = 3
        n = len(R)
        return {"ok": True, "chi2": chi2, "k": k, "n": n}
    except Exception:
        return {"ok": False}

def fit_mond(gal):
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon):
        Vbar2 = Vgas**2 + Upsilon * (Vdisk**2 + Vbulge**2)
        Vbar_m = np.sqrt(np.maximum(Vbar2, 0.0)) * 1000.0   # m/s
        R_m = R * 1e3 * pc
        gbar = np.maximum(Vbar_m**2 / R_m, 1e-40)
        y = gbar / a0_mond
        nu = 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-6))))
        return Vbar_m * np.sqrt(nu) / 1000.0

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=[0.5],
                            bounds=([0.1], [1.5]),
                            maxfev=30000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        k = 1
        n = len(R)
        return {"ok": True, "chi2": chi2, "k": k, "n": n}
    except Exception:
        return {"ok": False}

def fit_nfw_corrected(gal):
    """
    标准 NFW 拟合：使用 log(V200) 和 log(c) 作为内部参数。
    公式正确包含 c*x。
    """
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon, logV200, logc):
        V200 = np.exp(logV200)
        c = np.exp(logc)
        Vbar2 = Vgas**2 + Upsilon * (Vdisk**2 + Vbulge**2)
        R200 = V200 / (10.0 * H0_km_s_kpc)  # kpc
        x = np.maximum(R / R200, 1e-8)
        s = c * x                            # 关键修正
        f = np.log(1.0 + s) - s / (1.0 + s)  # 使用 s = c*x
        fc = np.log(1.0 + c) - c / (1.0 + c)
        Vhalo2 = V200**2 * f / (x * fc)
        return np.sqrt(np.maximum(Vbar2 + Vhalo2, 0.0))

    p0 = [0.5, np.log(120.0), np.log(10.0)]
    bounds = ([0.1, np.log(20.0), np.log(1.0)], [1.5, np.log(500.0), np.log(40.0)])

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=p0, bounds=bounds, maxfev=30000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        k = 3
        n = len(R)
        return {"ok": True, "chi2": chi2, "k": k, "n": n}
    except Exception:
        return {"ok": False}

# ----------------------------------------------------------------------
# 信息准则
# ----------------------------------------------------------------------
def aicc(chi2, n, k):
    if n <= k + 1:
        return np.inf
    return chi2 + 2*k + (2*k*(k+1)) / (n - k - 1)

def bic(chi2, n, k):
    return chi2 + k * np.log(n)

# ----------------------------------------------------------------------
# 主循环
# ----------------------------------------------------------------------
results = {"NEA": [], "MOND": [], "NFW": []}

for gal in galaxies:
    r_nea = fit_nea_smooth(gal)
    r_mond = fit_mond(gal)
    r_nfw = fit_nfw_corrected(gal)

    if not (r_nea["ok"] and r_mond["ok"] and r_nfw["ok"]):
        continue

    for name, r in [("NEA", r_nea), ("MOND", r_mond), ("NFW", r_nfw)]:
        n = r["n"]
        k = r["k"]
        chi2 = r["chi2"]
        r["redchi2"] = chi2 / (n - k)
        r["AICc"] = aicc(chi2, n, k)
        r["BIC"] = bic(chi2, n, k)
        results[name].append(r)

n_common = len(results["NEA"])
print(f"\nSuccessfully fitted all three models for {n_common} galaxies.")

# ----------------------------------------------------------------------
# 汇总
# ----------------------------------------------------------------------
print("\n" + "=" * 80)
print("模型表现对比（统一成功样本，修正NFW）")
print("=" * 80)

for name in ["NEA", "MOND", "NFW"]:
    arr = results[name]
    reds = np.array([r["redchi2"] for r in arr])
    aiccs = np.array([r["AICc"] for r in arr])
    bics = np.array([r["BIC"] for r in arr])
    print(f"\n{name}:")
    print(f"  median reduced chi2 = {np.median(reds):.3f}")
    print(f"  total AICc          = {np.sum(aiccs):.1f}")
    print(f"  total BIC           = {np.sum(bics):.1f}")
    print(f"  mean AICc           = {np.mean(aiccs):.1f}")
    print(f"  mean BIC            = {np.mean(bics):.1f}")

# 胜出比例（按 BIC 最低）
wins = []
for i in range(n_common):
    b = [results["NEA"][i]["BIC"], results["MOND"][i]["BIC"], results["NFW"][i]["BIC"]]
    wins.append(np.argmin(b))

wins = np.array(wins)
print("\n胜出比例（按 BIC 最低）:")
print(f"  NEA  : {np.sum(wins==0)} ({np.mean(wins==0)*100:.1f}%)")
print(f"  MOND : {np.sum(wins==1)} ({np.mean(wins==1)*100:.1f}%)")
print(f"  NFW  : {np.sum(wins==2)} ({np.mean(wins==2)*100:.1f}%)")