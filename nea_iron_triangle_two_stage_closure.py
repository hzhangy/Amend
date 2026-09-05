#!/usr/bin/env python3
# ==============================================================================
# N.E.A. Iron Triangle Two-Stage Closure Audit
# Strict version: compute all three raw values from CODATA + N.E.A. formulas
# ==============================================================================

import numpy as np

print("=" * 88)
print("  N.E.A. Iron Triangle Strict Two-Stage Closure Audit")
print("=" * 88)

# ------------------------------------------------------------------------------
# 1. CODATA constants
# ------------------------------------------------------------------------------
c = 299792458.0
G = 6.67430e-11
hbar = 1.054571817e-34

m_p = 1.67262192369e-27
m_e = 9.1093837015e-31

# ------------------------------------------------------------------------------
# 2. N.E.A. topological constants
# ------------------------------------------------------------------------------
Delta = 1.0 - np.sqrt(3.0) / 2.0
R = 1.0 / (1.0 + np.pi)
f_geo = 1.0 + Delta / (4.0 * np.pi)
N_max = np.exp(10.0 * np.sqrt(3.0))

# ------------------------------------------------------------------------------
# 3. Experimental dimensionless quantities
# ------------------------------------------------------------------------------
alphaG_inv_exp = hbar * c / (G * m_p**2)

# Z energy/mass anchor:
# Z_mass = m_e / (0.4*pi)
Z_mass = m_e / (0.4 * np.pi)

mpZ_exp = m_p / Z_mass

ell_P = np.sqrt(hbar * G / c**3)
lambda_Z = hbar / (Z_mass * c)
ellZ_exp = ell_P / lambda_Z

# Check exact physical identity
identity_exp = alphaG_inv_exp * ellZ_exp**2 * mpZ_exp**2

# ------------------------------------------------------------------------------
# 4. N.E.A. raw formulas
# ------------------------------------------------------------------------------
alphaG_inv_raw = N_max**5 / R
mpZ_raw = f_geo * np.sqrt(R) / (1.0 + R) * np.sqrt(N_max)
ellZ_raw = (1.0 + R) / (f_geo * N_max**3)

identity_raw = alphaG_inv_raw * ellZ_raw**2 * mpZ_raw**2

# ------------------------------------------------------------------------------
# 5. Two-stage correction mechanism
# ------------------------------------------------------------------------------
# Stage 1: volume-level correction
C1 = 1.0 + Delta**3 / 2.0

# Stage 2: solid-angle / null-mode correction
C2 = 1.0 + 3.0 * Delta**2 / (32.0 * np.pi**2)

# Allocation:
# alpha_G^-1 -> multiply by C1
# ell_P/lambda_Z -> multiply by C1^(-1/2) and C2
# m_p/Z -> multiply by C2^(-1)
alphaG_inv_c1 = alphaG_inv_raw * C1
ellZ_c1 = ellZ_raw * C1**(-0.5)
mpZ_c1 = mpZ_raw

alphaG_inv_final = alphaG_inv_c1
ellZ_final = ellZ_c1 * C2
mpZ_final = mpZ_c1 * C2**(-1.0)

identity_final = alphaG_inv_final * ellZ_final**2 * mpZ_final**2

# ------------------------------------------------------------------------------
# 6. Helper
# ------------------------------------------------------------------------------
def ppm(theory, exp):
    return (theory / exp - 1.0) * 1e6

def print_row(name, raw, c1, final, exp):
    print(
        f"{name:<20} | "
        f"{ppm(raw, exp):>16.3f} | "
        f"{ppm(c1, exp):>18.3f} | "
        f"{ppm(final, exp):>18.3f}"
    )

# ------------------------------------------------------------------------------
# 7. Report
# ------------------------------------------------------------------------------
print(f"[Delta]  = {Delta:.12f}")
print(f"[C1]     = 1 + Delta^3/2        = {C1:.12f}  ({(C1-1)*1e6:+.3f} ppm)")
print(f"[C2]     = 1 + 3 Delta^2/32pi^2 = {C2:.12f}  ({(C2-1)*1e6:+.3f} ppm)")
print("-" * 88)

print(f"{'Quantity':<20} | {'Raw ppm':>16} | {'After C1 ppm':>18} | {'After C1+C2 ppm':>18}")
print("-" * 88)

print_row("alpha_G^-1", alphaG_inv_raw, alphaG_inv_c1, alphaG_inv_final, alphaG_inv_exp)
print_row("ell_P/lambda_Z", ellZ_raw, ellZ_c1, ellZ_final, ellZ_exp)
print_row("m_p/Z", mpZ_raw, mpZ_c1, mpZ_final, mpZ_exp)

print("=" * 88)
print("Closure checks:")
print(f"  Experimental identity alphaG^-1 * ellZ^2 * mpZ^2 - 1 = {(identity_exp - 1)*1e6:+.6e} ppm")
print(f"  Raw N.E.A. identity                         - 1 = {(identity_raw - 1)*1e6:+.6e} ppm")
print(f"  Final corrected identity                    - 1 = {(identity_final - 1)*1e6:+.6e} ppm")
print("=" * 88)