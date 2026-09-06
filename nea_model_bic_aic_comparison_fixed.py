#!/usr/bin/env python3
"""
nea_model_bic_aic_comparison_fixed.py

Strict Model Selection on 143 SPARC Galaxies (Bug-Free NFW Implementation):
  1. N.E.A. Smooth Percolation (k=3: Upsilon, Rc, n)
  2. Standard MOND (k=1: Upsilon)
  3. Corrected Standard NFW Dark Matter Halo (k=3: Upsilon, V200, c)

Information Criteria:
  AICc = chi2 + 2k + 2k(k+1)/(n - k - 1)
  BIC  = chi2 + k * ln(n)
"""

import os
import glob
import numpy as np
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------
# Physical Constants
# ----------------------------------------------------------------------
pc = 3.08567758e16          # m
km = 1000.0                 # m
a0_mond = 1.2e-10           # m/s^2
H0_km_s_Mpc = 68.04
H0_km_s_kpc = H0_km_s_Mpc / 1000.0   # km/s / kpc

# ----------------------------------------------------------------------
# Data Loading (SPARC Database)
# ----------------------------------------------------------------------
files = sorted(glob.glob(os.path.join("Rotmod_LTG", "**", "*.dat"), recursive=True))
if not files:
    files = sorted(glob.glob("*.dat"))

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
    Vgas = d[:, 3]
    Vdisk = d[:, 4]
    Vbulge = d[:, 5] if d.shape[1] >= 6 else np.zeros_like(R)

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

print("=" * 80)
print(f"  Strict Model Selection on {len(galaxies)} SPARC Galaxies (Fixed NFW)")
print("=" * 80)

def get_eff_error(gal):
    return np.maximum(gal["dVobs"], 0.05 * gal["Vobs"])

# ----------------------------------------------------------------------
# 1. N.E.A. Smooth Model (k=3)
# ----------------------------------------------------------------------
def fit_nea_smooth(gal):
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon, Rc, n):
        Vbar2 = np.sign(Vgas)*(Vgas**2) + Upsilon * (Vdisk**2 + Vbulge**2)
        q = 1.0 + 1.0 / (1.0 + (R / Rc)**n)
        boost = (1.0 + R / Rc) ** (2.0 - q)
        return np.sqrt(np.maximum(Vbar2 * boost, 0.0))

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=[0.5, 5.0, 2.0],
                            bounds=([0.1, 0.5, 0.5], [1.5, 30.0, 6.0]),
                            maxfev=20000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        return {"ok": True, "chi2": chi2, "k": 3, "n": len(R), "params": popt}
    except Exception:
        return {"ok": False}

# ----------------------------------------------------------------------
# 2. Standard MOND (k=1)
# ----------------------------------------------------------------------
def fit_mond(gal):
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon):
        Vbar2 = np.sign(Vgas)*(Vgas**2) + Upsilon * (Vdisk**2 + Vbulge**2)
        Vbar_m = np.sqrt(np.maximum(Vbar2, 0.0)) * 1000.0
        R_m = R * 1e3 * pc
        gbar = np.maximum(Vbar_m**2 / R_m, 1e-40)
        y = gbar / a0_mond
        nu = 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-6))))
        return Vbar_m * np.sqrt(nu) / 1000.0

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=[0.5],
                            bounds=([0.1], [1.5]),
                            maxfev=20000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        return {"ok": True, "chi2": chi2, "k": 1, "n": len(R), "params": popt}
    except Exception:
        return {"ok": False}

# ----------------------------------------------------------------------
# 3. Corrected Standard NFW Halo Model (k=3)
# ----------------------------------------------------------------------
def fit_nfw_corrected(gal):
    R = gal["R"]
    Vobs = gal["Vobs"]
    dV = get_eff_error(gal)
    Vgas, Vdisk, Vbulge = gal["Vgas"], gal["Vdisk"], gal["Vbulge"]

    def func(R, Upsilon, V200, c):
        Vbar2 = np.sign(Vgas)*(Vgas**2) + Upsilon * (Vdisk**2 + Vbulge**2)
        R200 = V200 / (10.0 * H0_km_s_kpc)  # Virial radius in kpc
        x = np.maximum(R / R200, 1e-8)
        
        # CORRECT NFW DENSITY PROFILE INTEGRATION: s = r / r_s = c * x
        s = np.maximum(c * x, 1e-8)
        f = np.log(1.0 + s) - s / (1.0 + s)
        fc = np.log(1.0 + c) - c / (1.0 + c)
        
        Vhalo2 = (V200**2) * (1.0 / x) * (f / fc)
        return np.sqrt(np.maximum(Vbar2 + Vhalo2, 0.0))

    try:
        popt, _ = curve_fit(func, R, Vobs, sigma=dV,
                            p0=[0.5, 120.0, 8.0],
                            bounds=([0.1, 20.0, 1.0], [1.5, 400.0, 40.0]),
                            maxfev=20000)
        chi2 = np.sum(((Vobs - func(R, *popt)) / dV) ** 2)
        return {"ok": True, "chi2": chi2, "k": 3, "n": len(R), "params": popt}
    except Exception:
        return {"ok": False}

# ----------------------------------------------------------------------
# Information Criteria Functions
# ----------------------------------------------------------------------
def aicc(chi2, n, k):
    if n <= k + 1:
        return np.inf
    return chi2 + 2*k + (2*k*(k+1)) / (n - k - 1)

def bic(chi2, n, k):
    return chi2 + k * np.log(n)

# ----------------------------------------------------------------------
# Main Comparison Loop
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
print(f"Successfully fitted all three models for {n_common} galaxies.\n")

# ----------------------------------------------------------------------
# Output Final Results Table
# ----------------------------------------------------------------------
for name in ["NEA", "MOND", "NFW"]:
    arr = results[name]
    reds = np.array([r["redchi2"] for r in arr])
    aiccs = np.array([r["AICc"] for r in arr])
    bics = np.array([r["BIC"] for r in arr])
    print(f"{name}:")
    print(f"  median reduced chi2 = {np.median(reds):.3f}")
    print(f"  total AICc          = {np.sum(aiccs):.1f}")
    print(f"  total BIC           = {np.sum(bics):.1f}")
    print(f"  mean AICc           = {np.mean(aiccs):.1f}")
    print(f"  mean BIC            = {np.mean(bics):.1f}\n")

# Best Model by Lowest BIC
wins = []
for i in range(n_common):
    b = [results["NEA"][i]["BIC"], results["MOND"][i]["BIC"], results["NFW"][i]["BIC"]]
    wins.append(np.argmin(b))

wins = np.array(wins)
print("=" * 80)
print("Model Selection Win Rate (by Lowest BIC):")
print(f"  N.E.A. : {np.sum(wins==0)}/{n_common} ({np.mean(wins==0)*100:.1f}%)")
print(f"  MOND   : {np.sum(wins==1)}/{n_common} ({np.mean(wins==1)*100:.1f}%)")
print(f"  NFW    : {np.sum(wins==2)}/{n_common} ({np.mean(wins==2)*100:.1f}%)")
print("=" * 80)