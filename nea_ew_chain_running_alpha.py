"""
nea_ew_chain_running_alpha.py

Correct the electroweak chain by using the running alpha at M_Z
instead of the low-energy Thomson limit.

This isolates how much of the m_W / m_Z residual comes from the
wrong scale choice, versus the remaining v_h defect.
"""

import numpy as np

pi = np.pi
sqrt3 = np.sqrt(3.0)

R = 1.0 / (1.0 + pi)
Delta = 1.0 - sqrt3 / 2.0
delta = Delta / (4.0 * pi)

# N.E.A. alpha values
alpha_inv_0 = 25.0 * sqrt3 * pi + 1.0      # 137.035 (Thomson / IR)
alpha_inv_MZ_nea = 126.56                  # one-loop running with NEA masses
alpha_inv_MZ_exp = 127.9                   # PDG effective value at M_Z

# v_h from N.E.A. tree formula (unchanged)
m_p_GeV = 0.93827208816
v_h = m_p_GeV * alpha_inv_0 * 6.0 / pi

# sin^2 theta_W corrected by C4
C4 = Delta / (32.0 * pi**2)
s2w_raw = R - delta
s2w = s2w_raw + C4

def mW_from_alpha_inv(alpha_inv, s2w_val, vh_val):
    alpha = 1.0 / alpha_inv
    return vh_val * np.sqrt(pi * alpha / s2w_val)

def mZ_from_mW(mW_val, s2w_val):
    return mW_val / np.sqrt(1.0 - s2w_val)

mW_exp = 80.377
mW_err = 0.012
mZ_exp = 91.1876
mZ_err = 0.0021

print("v_h (fixed)   =", round(v_h, 4), "GeV")
print("sin^2 theta_W =", round(s2w, 6))
print()
print(f"{'alpha_inv used':<18} {'m_W (GeV)':>10} {'m_W sigma':>10} {'m_Z (GeV)':>10} {'m_Z sigma':>10}")
print("-" * 58)

for name, inv in [("alpha(0)=137.035", alpha_inv_0),
                  ("alpha(MZ)=127.9", alpha_inv_MZ_exp),
                  ("alpha(MZ)=126.56", alpha_inv_MZ_nea)]:
    mW = mW_from_alpha_inv(inv, s2w, v_h)
    mZ = mZ_from_mW(mW, s2w)
    pullW = (mW - mW_exp) / mW_err
    pullZ = (mZ - mZ_exp) / mZ_err
    print(f"{name:<18} {mW:>10.3f} {pullW:>+10.1f} {mZ:>10.3f} {pullZ:>+10.1f}")