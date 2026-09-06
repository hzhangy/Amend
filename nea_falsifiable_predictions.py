#!/usr/bin/env python3
"""
nea_falsifiable_predictions.py

严格验证 N.E.A. 的两个正向可证伪实验：
  1. 高红移临界表面密度的宇宙学演化
     Sigma_crit(z) = Sigma_crit(0) * H(z)/H_0
  2. UDG 外场过渡宽度的数量级展宽
     delta_iso / delta_cluster ~ Sigma_cluster / Sigma_isolated

用 CODATA 2018 常数与 Planck 2018 宇宙学参数严格计算。
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# ============================================================
# 常数
# ============================================================
c       = 2.99792458e8       # m/s
G       = 6.67430e-11        # m^3 kg^-1 s^-2
H0_si   = 2.20495e-18        # s^-1  (68.04 km/s/Mpc)
H0_planck = 2.1758e-18       # s^-1  (67.4 km/s/Mpc)

# 宇宙学参数 (Planck 2018)
Omega_m = 0.3153
Omega_L = 0.6847

# N.E.A. 拓扑常数
pi    = np.pi
R     = 1.0 / (1.0 + pi)
Delta = 1.0 - np.sqrt(3.0) / 2.0

# 单位转换
M_sun   = 1.98892e30    # kg
pc      = 3.08567758e16 # m
kg_per_m2_per_Msun_pc2 = M_sun / pc**2

print("=" * 78)
print("  N.E.A. 正向可证伪实验严格验证")
print("=" * 78)

# ============================================================
# 1. 临界表面密度裸值
# ============================================================
Sigma_crit_0 = c * H0_si / (8.0 * pi * G)  # kg/m^2
Sigma_crit_0_solar = Sigma_crit_0 / kg_per_m2_per_Msun_pc2

print(f"\n【1】裸临界表面密度 Σ_crit(0)")
print(f"  {Sigma_crit_0:.6f} kg/m^2 = {Sigma_crit_0_solar:.2f} M_sun/pc^2")

# ============================================================
# 2. (1+R) 理论候选修正
# ============================================================
dressed_theory = Sigma_crit_0_solar * (1.0 + R)
dressed_empir  = Sigma_crit_0_solar * 1.350

print(f"\n【2】理论 vs 经验 dressed 因子")
print(f"  (1+R) 理论候选: {dressed_theory:.2f} M_sun/pc^2")
print(f"  SPARC 经验拟合: {dressed_empir:.2f} M_sun/pc^2")
print(f"  偏差: {abs(dressed_empir-dressed_theory)/dressed_theory*100:.1f}%")

# ============================================================
# 3. 星系团有效 q 值
# ============================================================
def q_logistic(Sigma, Sigma_c):
    """logistic 形式 q(Σ)，不是 G 篇原公式，仅作数值估计"""
    return 1.0 + np.tanh(Sigma / Sigma_c)

Sigma_coma = 141.0  # M_sun/pc^2, Coma 团平均面密度
x_coma = Sigma_coma / Sigma_crit_0_solar
q_coma = q_logistic(Sigma_coma, Sigma_crit_0_solar)

Sigma_cluster_center = 700.0  # M_sun/pc^2, 星系团中心
q_center = q_logistic(Sigma_cluster_center, Sigma_crit_0_solar)

Sigma_cluster_outer = 50.0   # M_sun/pc^2, 星系团外围
q_outer = q_logistic(Sigma_cluster_outer, Sigma_crit_0_solar)

print(f"\n【3】星系团有效维度 q")
print(f"  中心 (Σ~700): q = {q_center:.3f}")
print(f"  平均 (Σ~141): q = {q_coma:.3f}")
print(f"  外围 (Σ~50):  q = {q_outer:.3f}")

# ============================================================
# 4. 高红移 Σ_crit(z) = Σ_crit(0) * H(z)/H_0
# ============================================================
def H_ratio(z):
    """H(z)/H0 for flat LambdaCDM"""
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def Sigma_crit_z(z):
    """Σ_crit(z) in M_sun/pc^2"""
    return Sigma_crit_0_solar * H_ratio(z)

print(f"\n【4】高红移临界表面密度演化")
print(f"  {'z':>5} {'H(z)/H0':>10} {'Σ_crit(z) [M_sun/pc^2]':>25}")
z_test = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
for z in z_test:
    hr = H_ratio(z)
    sc = Sigma_crit_z(z)
    print(f"  {z:5.1f} {hr:10.4f} {sc:25.1f}")

# ============================================================
# 5. 正向可证伪实验 I 的判决标准
# ============================================================
# 设一个星系面密度为 300 M_sun/pc^2
Sigma_gal = 300.0  # M_sun/pc^2

print(f"\n【5】正向可证伪实验 I：星系面密度 300 M_sun/pc^2")
print(f"  在 z=0: Σ/Σ_crit = {Sigma_gal/Sigma_crit_0_solar:.2f} → 牛顿区 (q≈2)")
for z in [1.0, 2.0, 3.0]:
    sc = Sigma_crit_z(z)
    ratio = Sigma_gal / sc
    q_gal = q_logistic(Sigma_gal, sc)
    print(f"  在 z={z}: Σ/Σ_crit = {ratio:.2f} → q = {q_gal:.3f}")

# ============================================================
# 6. 正向可证伪实验 II：UDG 过渡宽度比
# ============================================================
# δ_0 ∝ 1/Σ_background 的简单模型
Sigma_isolated = 10.0    # M_sun/pc^2，极端孤立 UDG
Sigma_cluster  = 100.0   # M_sun/pc^2，星系团背景

delta_ratio = Sigma_cluster / Sigma_isolated

print(f"\n【6】正向可证伪实验 II：UDG 过渡宽度比")
print(f"  Σ_isolated = {Sigma_isolated} M_sun/pc^2")
print(f"  Σ_cluster  = {Sigma_cluster} M_sun/pc^2")
print(f"  δ_iso/δ_cluster ~ Σ_cluster/Σ_isolated = {delta_ratio:.1f}")
print(f"  预言：孤立 UDG 过渡带宽是星系团内 UDG 的约 {delta_ratio:.0f} 倍")

# ============================================================
# 7. 绘图：高红移临界密度演化 + 判决区
# ============================================================
z_arr = np.linspace(0, 3.5, 200)
sc_z = np.array([Sigma_crit_z(z) for z in z_arr])
hr_arr = np.array([H_ratio(z) for z in z_arr])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=130)

# 左图：H(z)/H0
ax1.plot(z_arr, hr_arr, 'b-', linewidth=2.2, label='$H(z)/H_0$')
ax1.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='$z=0$ reference')
ax1.set_xlabel('Redshift $z$', fontsize=12)
ax1.set_ylabel('$H(z)/H_0$', fontsize=12)
ax1.set_title('Cosmological Expansion Rate', fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)

# 右图：Σ_crit(z)
ax2.plot(z_arr, sc_z, 'r-', linewidth=2.5, label='$\\Sigma_{\\rm crit}(z)$')
ax2.axhline(y=Sigma_crit_0_solar, color='k', linestyle='--', alpha=0.5, label='$\\Sigma_{\\rm crit}(0)$')
ax2.axhline(y=Sigma_gal, color='g', linestyle=':', linewidth=1.5, label='Galaxy $\\Sigma=300\\ M_\\odot/\\rm pc^2$')

# 标注判决点
for z_test_pt in [1.0, 2.0]:
    sc_pt = Sigma_crit_z(z_test_pt)
    ax2.plot(z_test_pt, sc_pt, 'o', markersize=8, color='purple')

ax2.set_xlabel('Redshift $z$', fontsize=12)
ax2.set_ylabel('$\\Sigma_{\\rm crit}(z)\\ [M_\\odot/\\rm pc^2]$', fontsize=12)
ax2.set_title('High-Redshift Critical Surface Density', fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)

plt.tight_layout()
plt.savefig('nea_falsifiable_predictions.png', dpi=150)
print("\n>> [完成] 图表已保存为: nea_falsifiable_predictions.png")
plt.show()

print("=" * 78)
print("  验证完成")
print("=" * 78)