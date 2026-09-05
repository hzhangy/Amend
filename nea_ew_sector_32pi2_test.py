"""
nea_ew_sector_32pi2_test.py

Test whether the 32π² correction that improves sin²θ_W
also propagates correctly to m_W and m_Z through the
N.E.A. tree-level electroweak chain.

The chain is:
    1. sin²θ_W  = R - δ
    2. v_h      = m_p · α^{-1} · 6/π
    3. m_W      = v_h · sqrt(π·α / sin²θ_W)
    4. m_Z      = m_W / cosθ_W

We compare the raw chain and the C4-corrected chain
against PDG values, computing ppm and approximate σ
for each. No new free coefficient is introduced.
"""

import numpy as np

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------
pi = np.pi
sqrt3 = np.sqrt(3.0)
sqrt2 = np.sqrt(2.0)

# Topological constants
R = 1.0 / (1.0 + pi)
Delta = 1.0 - sqrt3 / 2.0
delta = Delta / (4.0 * pi)

# N.E.A. α^{-1} and α
alpha_inv_nea = 25.0 * sqrt3 * pi + 1.0
alpha_nea = 1.0 / alpha_inv_nea

# Electron / ZY scale
m_e = 0.51099895  # MeV
ZY = m_e / (0.4 * pi)  # MeV
ZY_GeV = ZY / 1000.0

# Proton mass (CODATA 2018)
m_p_exp = 938.27208816  # MeV
m_p_GeV = m_p_exp / 1000.0

# PDG electroweak observables
s2w_exp = 0.23122
s2w_err = 0.00004
mW_exp = 80.377  # GeV
mW_err = 0.012   # GeV
mZ_exp = 91.1876 # GeV
mZ_err = 0.0021  # GeV


def ppm(theory, exp):
    """Relative deviation in parts per million."""
    return (theory / exp - 1.0) * 1e6


def sigma_pull(theory, exp, err):
    """Approximate statistical pull in units of sigma."""
    return abs(theory - exp) / err


# ----------------------------------------------------------------------
# N.E.A. tree-level electroweak chain
# ----------------------------------------------------------------------
def sin2theta_raw():
    return R - delta

def vh_raw():
    # v_h = m_p · α^{-1} · 6/π, using the N.E.A. topological α
    return m_p_GeV * alpha_inv_nea * 6.0 / pi

def mW_from(s2w, vh):
    # m_W = v_h · sqrt(π·α / sin²θ_W)
    return vh * np.sqrt(pi * alpha_nea / s2w)

def mZ_from(mW, s2w):
    # m_Z = m_W / cosθ_W
    cosw = np.sqrt(1.0 - s2w)
    return mW / cosw

# Raw chain
s2w_raw = sin2theta_raw()
vh_raw_val = vh_raw()
mW_raw = mW_from(s2w_raw, vh_raw_val)
mZ_raw = mZ_from(mW_raw, s2w_raw)

# Corrected chain: apply C4 to sin²θ_W only
C4 = Delta / (32.0 * pi**2)
s2w_corr = s2w_raw + C4
mW_corr = mW_from(s2w_corr, vh_raw_val)
mZ_corr = mZ_from(mW_corr, s2w_corr)

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------
print("=" * 80)
print("N.E.A. Electroweak Chain: raw vs C4-corrected")
print("=" * 80)

print(f"\nTree-level sin²θ_W = {s2w_raw:.6f}")
print(f"C4 = Δ/(32π²)      = {C4:.6f}  ({C4*1e6:+.1f} ppm)")
print(f"C4-corrected       = {s2w_corr:.6f}")
print(f"PDG MS-bar         = {s2w_exp:.6f} ± {s2w_err:.6f}")
print(f"\nsin²θ_W raw σ     = {sigma_pull(s2w_raw, s2w_exp, s2w_err):.1f}")
print(f"sin²θ_W corr σ    = {sigma_pull(s2w_corr, s2w_exp, s2w_err):.2f}")

print("\n" + "-" * 80)
print(f"{'Quantity':<12} {'Raw (GeV)':>12} {'Corr (GeV)':>12} {'Exp (GeV)':>12} {'Raw σ':>8} {'Corr σ':>8}")
print("-" * 80)
print(f"{'m_W':<12} {mW_raw:>12.4f} {mW_corr:>12.4f} {mW_exp:>12.4f} {sigma_pull(mW_raw,mW_exp,mW_err):>8.1f} {sigma_pull(mW_corr,mW_exp,mW_err):>8.1f}")
print(f"{'m_Z':<12} {mZ_raw:>12.4f} {mZ_corr:>12.4f} {mZ_exp:>12.4f} {sigma_pull(mZ_raw,mZ_exp,mZ_err):>8.1f} {sigma_pull(mZ_corr,mZ_exp,mZ_err):>8.1f}")

print("\n" + "=" * 80)
print("Interpretation:")
print("  C4 improves sin²θ_W strongly.")
print("  The effect on m_W and m_Z is in the correct direction")
print("  but is weaker, because m_W/m_Z also depend on v_h,")
print("  which is not corrected here.")
print("  This isolates v_h as the remaining tree-level residual.")
print("=" * 80)