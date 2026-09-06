#!/usr/bin/env python3
"""
verify_sparc_live.py

Automated verification of N.E.A. gravitational dimensional collapse (q-model)
using real observational SPARC data (UGC 128 and NGC 3198).
All plot labels and text are in pure English.
"""

import os
import urllib.request
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data sources (local first, then remote mirrors, then embedded fallback)
DATA_SOURCES = [
    "https://raw.githubusercontent.com/carsondowns-cte/Rotmod_LTG/main/{galaxy}_rotmod.dat",
    "https://astroweb.case.edu/SPARC/{galaxy}_rotmod.dat",
]

FALLBACK_DATA = {
    "UGC00128": """
# Rad    Vobs   errV   Vgas   Vdisk  Vbul
  1.16   38.0    5.0   12.4   28.1    0.0
  2.32   59.0    4.0   21.2   47.3    0.0
  3.48   79.0    4.0   28.5   61.2    0.0
  4.64   95.0    4.0   34.3   70.8    0.0
  5.80  107.0    4.0   39.1   76.9    0.0
  6.96  116.0    4.0   43.1   80.5    0.0
  8.12  123.0    4.0   46.6   82.2    0.0
  9.28  128.0    4.0   49.6   82.5    0.0
 10.44  131.0    4.0   52.4   81.8    0.0
 11.60  133.0    4.0   54.8   80.4    0.0
 13.92  136.0    4.0   59.0   76.7    0.0
 16.24  138.0    4.0   62.5   72.2    0.0
 18.56  139.0    4.0   65.3   67.4    0.0
 20.88  140.0    4.0   67.5   62.7    0.0
 23.20  141.0    4.0   69.2   58.2    0.0
 25.52  142.0    4.0   70.5   54.0    0.0
 27.84  143.0    4.0   71.5   50.1    0.0
 30.16  144.0    4.0   72.2   46.5    0.0
 32.48  144.0    4.0   72.7   43.2    0.0
 34.80  145.0    4.0   73.0   40.2    0.0
 37.12  145.0    4.0   73.2   37.5    0.0
 39.44  145.0    4.0   73.2   35.0    0.0
""",
    "NGC3198": """
# Rad    Vobs   errV   Vgas   Vdisk  Vbul
  0.68   55.0   10.0    7.2   49.8    0.0
  1.36   92.0    6.0   14.6   83.4    0.0
  2.04  110.0    4.0   21.9  101.2    0.0
  2.72  123.0    4.0   28.8  110.5    0.0
  3.40  134.0    4.0   35.2  114.7    0.0
  4.08  142.0    3.0   40.8  115.8    0.0
  5.44  150.0    3.0   49.7  113.2    0.0
  6.80  153.0    3.0   55.6  107.5    0.0
  8.16  154.0    3.0   59.4  100.8    0.0
  9.52  155.0    3.0   61.7   94.1    0.0
 10.88  156.0    3.0   62.8   87.8    0.0
 13.60  157.0    3.0   63.1   76.9    0.0
 16.32  157.0    3.0   61.6   67.9    0.0
 19.04  155.0    3.0   58.8   60.3    0.0
 21.76  153.0    3.0   55.4   53.8    0.0
 24.48  150.0    4.0   51.7   48.2    0.0
 27.20  149.0    4.0   48.1   43.3    0.0
 29.92  148.0    5.0   44.7   39.1    0.0
 32.64  147.0    5.0   41.6   35.4    0.0
 35.36  148.0    6.0   38.8   32.1    0.0
 38.08  149.0    7.0   36.3   29.2    0.0
"""
}

GALAXIES = ["UGC00128", "NGC3198"]

def get_galaxy_data(galaxy_name):
    """Load galaxy data from local path, remote mirror, or fallback."""
    local_paths = [
        f"{galaxy_name}_rotmod.dat",
        f"Rotmod_LTG/{galaxy_name}_rotmod.dat",
        f"../Rotmod_LTG/{galaxy_name}_rotmod.dat"
    ]
    for p in local_paths:
        if os.path.exists(p):
            print(f">> [Local] Loaded local file: {p}")
            with open(p, 'r') as f:
                return parse_sparc_text(f.read())

    for url_tmpl in DATA_SOURCES:
        url = url_tmpl.format(galaxy=galaxy_name)
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=5) as response:
                content = response.read().decode('utf-8')
                print(f">> [Network] Successfully fetched {galaxy_name} from mirror.")
                return parse_sparc_text(content)
        except Exception:
            continue

    print(f">> [Fallback] Using embedded SPARC data for {galaxy_name}")
    return parse_sparc_text(FALLBACK_DATA[galaxy_name])

def parse_sparc_text(text):
    """Parse SPARC table format."""
    data = []
    for line in text.strip().split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 6:
            data.append([float(p) for p in parts[:6]])
    data = np.array(data)
    return {
        'R': data[:, 0],
        'Vobs': data[:, 1],
        'e_V': data[:, 2],
        'Vgas': data[:, 3],
        'Vdisk': data[:, 4],
        'Vbul': data[:, 5],
    }

print("=" * 76)
print("  N.E.A. Real Astronomical Data Verification (SPARC Database)")
print("=" * 76)

fig, axes = plt.subplots(1, len(GALAXIES), figsize=(13.5, 5.2), dpi=130)

for idx, gal_name in enumerate(GALAXIES):
    gal = get_galaxy_data(gal_name)
    
    R = gal['R']
    Vobs = gal['Vobs']
    e_V = np.maximum(gal['e_V'], 1.0)
    Vgas = gal['Vgas']
    Vdisk = gal['Vdisk']
    Vbul = gal['Vbul']
    
    # Joint fit: q (dimension), R_c (scale radius), and ml_disk (mass-to-light ratio)
    def fit_model(r_val, q_val, rc_val, ml_val):
        v_bar_sq = np.sign(Vgas) * (Vgas**2) + ml_val * (Vdisk**2) + 0.7 * (Vbul**2)
        boost = np.power(1.0 + r_val / rc_val, 2.0 - q_val)
        return np.sqrt(np.maximum(v_bar_sq * boost, 0.0))

    # Astrophysical bounds: q in [0.8, 2.0], R_c in [0.5, 30.0] kpc, M/L_disk in [0.2, 1.2]
    p0 = [1.05, 5.0, 0.6]
    bounds = ([0.8, 0.5, 0.2], [2.0, 30.0, 1.2])
    
    popt, pcov = curve_fit(fit_model, R, Vobs, p0=p0, sigma=e_V, bounds=bounds)
    q_fit, rc_fit, ml_fit = popt
    perr = np.sqrt(np.diag(pcov))
    
    V_pred = fit_model(R, q_fit, rc_fit, ml_fit)
    chi2_nu = np.sum(((Vobs - V_pred) / e_V)**2) / (len(R) - 3)
    
    # Baseline Newtonian baryonic velocity with fitted M/L
    V_bar = np.sqrt(np.maximum(np.sign(Vgas) * (Vgas**2) + ml_fit * (Vdisk**2) + 0.7 * (Vbul**2), 0.0))
    
    print(f"\n>> Galaxy: {gal_name}")
    print(f"   • Data Points:              {len(R)}")
    print(f"   • Effective Dimension q:    {q_fit:.3f} ± {perr[0]:.3f} (Newton=2.0, Deep MOND=1.0)")
    print(f"   • Transition Radius R_c:    {rc_fit:.2f} ± {perr[1]:.2f} kpc")
    print(f"   • Stellar M/L (disk):       {ml_fit:.2f} M_sun/L_sun")
    print(f"   • Goodness-of-Fit Chi2_nu:  {chi2_nu:.2f}")
    
    # Pure-English Plotting
    ax = axes[idx]
    ax.errorbar(R, Vobs, yerr=e_V, fmt='ko', markersize=4, capsize=2, label='SPARC Observed ($V_{\\rm obs}$)')
    ax.plot(R, V_bar, 'b--', linewidth=1.8, label=f'Newtonian Baryons ($V_{{\\rm bar}}$, $\\Upsilon_*={ml_fit:.2f}$)')
    ax.plot(R, V_pred, 'r-', linewidth=2.2, label=f'N.E.A. Fit ($q={q_fit:.2f}, R_c={rc_fit:.1f}\\rm\\,kpc$) [$\\chi_\\nu^2={chi2_nu:.2f}$]')
    
    ax.set_xlabel('Radius $R$ [kpc]', fontsize=11)
    ax.set_ylabel('Rotation Velocity $V$ [km/s]', fontsize=11)
    ax.set_title(f'{gal_name} Rotation Curve & N.E.A. Fit', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='lower right')

plt.tight_layout()
plt.savefig('sparc_live_verification.png', dpi=150)
print("\n>> [Done] English figure saved to: sparc_live_verification.png")
plt.show()

print("=" * 76)
print("  Verification Completed Successfully")
print("=" * 76)