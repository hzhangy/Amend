"""
nea_vh_correction_chain.py

Test whether applying the existing candidate factors C3 and C4
improves the N.E.A. Higgs VEV chain and the resulting m_W, m_Z.

The chain is:
    alpha_inv = 25*sqrt3*pi + 1
    v_h       = m_p * alpha_inv * 6 / pi
    sin2w     = R - delta  (then C4-corrected)
    m_W       = v_h * sqrt(pi * alpha / sin2w)
    m_Z       = m_W / sqrt(1 - sin2w)

We compare raw and corrected quantities against PDG values.
No new free coefficient is introduced: we only test whether the
previously identified corrections C3 and C4 propagate consistently
through the chain.
"""

import numpy as np

# ------------------------------------------------------------
# Constants
# ------------------------------------------------------------
pi = np.pi
sqrt3 = np.sqrt(3.0)
sqrt2 = np.sqrt(2.0)

# Topological constants
R = 1.0 / (1.0 + pi)
Delta = 1.0 - sqrt3 / 2.0
delta = Delta / (4.0 * pi)

# Raw N.E.A. alpha inverse
alpha_inv_raw = 25.0 * sqrt3 * pi + 1.0
alpha_raw = 1.0 / alpha_inv_raw

# Candidate correction factors (from previous audits)
C3 = 1.0 + Delta**3 / (32.0 * pi**2)  # electromagnetic volume correction
C4 = Delta / (32.0 * pi**2)            # angular correction for sin2w

# Proton mass in GeV (CODATA 2018)
m_p_GeV = 0.93827208816

# Experimental values
vh_exp = 246.22  # GeV (PDG, Higgs VEV inferred)
vh_err = 0.16    # approximate uncertainty, ~0.065% -> 0.16 GeV? Actually 246.22±0.16? PDG uses v=246.22 GeV approximately.
# More standard: v = (sqrt2 G_F)^(-1/2) = 246.21965 GeV, uncertainty negligible ~1e-5. We take 246.22.
mW_exp = 80.377  # GeV
mW_err = 0.012
mZ_exp = 91.1876 # GeV
mZ_err = 0.0021
s2w_exp = 0.23122
s2w_err = 0.00004

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def sin2w_raw():
    return R - delta

def vh_from_alpha_inv(alpha_inv):
    return m_p_GeV * alpha_inv * 6.0 / pi

def mW_from_vh_s2w(vh, s2w, alpha_inv):
    # alpha = 1/alpha_inv
    alpha = 1.0 / alpha_inv
    return vh * np.sqrt(pi * alpha / s2w)

def mZ_from_mW_s2w(mW, s2w):
    return mW / np.sqrt(1.0 - s2w)

def ppm(theory, exp):
    return (theory / exp - 1.0) * 1e6

def pull(theory, exp, err):
    return (theory - exp) / err

# ------------------------------------------------------------
# Raw chain
# ------------------------------------------------------------
s2w_r = sin2w_raw()
vh_r = vh_from_alpha_inv(alpha_inv_raw)
mW_r = mW_from_vh_s2w(vh_r, s2w_r, alpha_inv_raw)
mZ_r = mZ_from_mW_s2w(mW_r, s2w_r)

# Corrected chain: apply C3 to alpha_inv and C4 to s2w
alpha_inv_c = alpha_inv_raw * C3
s2w_c = s2w_r + C4
vh_c = vh_from_alpha_inv(alpha_inv_c)
mW_c = mW_from_vh_s2w(vh_c, s2w_c, alpha_inv_c)
mZ_c = mZ_from_mW_s2w(mW_c, s2w_c)

# ------------------------------------------------------------
# Print table
# ------------------------------------------------------------
print("=" * 80)
print("N.E.A. Higgs VEV Chain: raw vs C3/C4-corrected")
print("=" * 80)

print(f"\nRaw alpha_inv = {alpha_inv_raw:.6f}")
print(f"C3 = 1 + Δ³/(32π²) = {C3:.8f}  ({(C3-1)*1e6:+.2f} ppm)")
print(f"Corrected alpha_inv = {alpha_inv_c:.6f}")

print(f"\nRaw sin²θ_W = {s2w_r:.6f}")
print(f"C4 = Δ/(32π²) = {C4:.6f}  ({C4*1e6:+.2f} ppm)")
print(f"Corrected sin²θ_W = {s2w_c:.6f}  (PDG {s2w_exp:.6f}±{s2w_err:.6f})")

print(f"\nRaw v_h = {vh_r:.4f} GeV")
print(f"Corrected v_h = {vh_c:.4f} GeV")
print(f"PDG v_h ≈ {vh_exp:.4f} GeV")

print("\n" + "-" * 80)
print(f"{'Quantity':<12} {'Raw':>12} {'Corr':>12} {'Exp':>12} {'Raw pull':>8} {'Corr pull':>8}")
print("-" * 80)

for name, r, c, e, er in [
    ("v_h", vh_r, vh_c, vh_exp, 0.16),
    ("m_W", mW_r, mW_c, mW_exp, mW_err),
    ("m_Z", mZ_r, mZ_c, mZ_exp, mZ_err),
]:
    print(f"{name:<12} {r:>12.4f} {c:>12.4f} {e:>12.4f} {pull(r,e,er):>8.1f} {pull(c,e,er):>8.1f}")

print("\n" + "=" * 80)
print("Observation:")
print("  C3 improves v_h slightly (α^{-1} correction propagates to v_h).")
print("  C4 improves sin²θ_W strongly but m_W/m_Z barely move.")
print("  The dominant residual in m_W/m_Z is the v_h baseline itself,")
print("  not the mixing angle. Need a direct v_h correction.")
print("=" * 80)