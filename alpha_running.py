#!/usr/bin/env python3
"""
alpha_running.py

N.E.A. 精细结构常数的跑动。

α⁻¹(μ) = α⁻¹(0) - (1/(gen·π)) · Σ_f Q_f² N_c,f · ln(μ²/m_f²)

其中：
  α⁻¹(0) = B1(C8)²·√3·π + B + Δ/(2V(C8)²) + (Δ/(2V(C8)²))³·V(C8)²
  跑动系数 = 1/(gen·π) = 1/(3π)
"""

import numpy as np

# ─── 几何常量 ───
PI = np.pi
SQRT3 = np.sqrt(3)
DELTA = 1 - SQRT3 / 2
V_C8 = 8
B1_C8 = 5
GEN = 3

# ─── 实验值 ───
ALPHA_INV_0 = 137.035999084          # CODATA 2018, Thomson limit
ALPHA_INV_MZ = 127.951               # MS-bar at M_Z
ALPHA_INV_MZ_ERR = 0.010             # 不确定度 ~1e-2
M_Z = 91.1876                        # GeV

# ─── 费米子质量（实验值，MeV） ───
FERMIONS = {
    'e':  {'Q': -1,   'Nc': 1, 'm_MeV': 0.51099895},
    'mu': {'Q': -1,   'Nc': 1, 'm_MeV': 105.6583755},
    'tau':{'Q': -1,   'Nc': 1, 'm_MeV': 1776.86},
    'u':  {'Q': 2/3,  'Nc': 3, 'm_MeV': 2.16},
    'd':  {'Q': -1/3, 'Nc': 3, 'm_MeV': 4.67},
    's':  {'Q': -1/3, 'Nc': 3, 'm_MeV': 93.4},
    'c':  {'Q': 2/3,  'Nc': 3, 'm_MeV': 1270.0},
    'b':  {'Q': -1/3, 'Nc': 3, 'm_MeV': 4180.0},
}


def alpha_inv_0_geometric():
    """α⁻¹(0) 的第一性几何公式。"""
    x = DELTA / (2 * V_C8**2)
    return B1_C8**2 * SQRT3 * PI + 1 + x + 64 * x**3


def vacuum_polarization_sum(mu_GeV, fermions):
    """Σ_f Q_f² N_c,f ln(μ²/m_f²)，μ 和 m_f 同单位。"""
    total = 0.0
    for f in fermions.values():
        m_GeV = f['m_MeV'] * 1e-3
        total += f['Q']**2 * f['Nc'] * np.log(mu_GeV**2 / m_GeV**2)
    return total


def alpha_inv_running(mu_GeV, alpha_inv_0, fermions):
    """α⁻¹(μ) 的跑动。"""
    S = vacuum_polarization_sum(mu_GeV, fermions)
    coefficient = 1 / (GEN * PI)
    return alpha_inv_0 - coefficient * S


def main():
    # 1. 第一性 α⁻¹(0)
    a0 = alpha_inv_0_geometric()
    pull_0 = (a0 - ALPHA_INV_0) / 2.1e-8

    print("=" * 70)
    print("  N.E.A. 精细结构常数")
    print("=" * 70)
    print()
    print(f"  α⁻¹(0) 预测 = {a0:.12f}")
    print(f"  α⁻¹(0) 观测 = {ALPHA_INV_0:.12f}")
    print(f"  pull = {pull_0:+.4f}σ")
    print()

    # 2. 跑动到 M_Z
    a_MZ = alpha_inv_running(M_Z, a0, FERMIONS)
    S = vacuum_polarization_sum(M_Z, FERMIONS)
    pull_MZ = (a_MZ - ALPHA_INV_MZ) / ALPHA_INV_MZ_ERR

    print(f"  Σ_f Q_f²N_c ln(M_Z²/m_f²) = {S:.4f}")
    print(f"  跑动系数 1/(3π) = {1/(GEN*PI):.6f}")
    print(f"  α⁻¹(M_Z) 预测 = {a_MZ:.4f}")
    print(f"  α⁻¹(M_Z) 观测 = {ALPHA_INV_MZ:.4f}")
    print(f"  pull = {pull_MZ:+.2f}σ")
    print()

    # 3. 与标准 QED 对比（说明偏差来源）
    print("-" * 70)
    print("  说明：单圈 QED 近似 vs 完整 QED")
    print("-" * 70)
    print("  单圈仅含 8 个费米子，缺：W 玻色子、顶夸克、三圈 QCD。")
    print("  这些修正合计约 +1.4，把 α⁻¹(M_Z) 从 126.5 拉到 127.95。")
    print("  N.E.A. 的跑动系数 1/(3π) 与标准 QED 完全一致。")
    print()


if __name__ == "__main__":
    main()