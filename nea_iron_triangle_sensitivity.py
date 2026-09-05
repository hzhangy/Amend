"""
nea_iron_triangle_sensitivity.py

解析验证：尺度骨架铁三角的对数全微分敏感度矩阵与残差联合消除分析
"""

import numpy as np

# 1. 基础拓扑参数定义
U_w = 10.0 * np.sqrt(3.0)  # 弱激活租金 ≈ 17.320508
R = 1.0 / (1.0 + np.pi)  # 投影残差 ≈ 0.241523
Delta = 1.0 - np.sqrt(3.0) / 2.0  # 四面体间隙 ≈ 0.133975
f_g = 1.0 + Delta / (4.0 * np.pi)  # 一阶立体角修正 ≈ 1.010660

# 2. 实验基准目标值 (CODATA 2018)
# G = 6.67430e-11, m_p = 1.67262192369e-27, hbar = 1.054571817e-34, c = 299792458
Y1_target = (1.054571817e-34 * 299792458) / (
    6.67430e-11 * (1.67262192369e-27) ** 2
)  # alpha_G^-1 ≈ 1.693054e38
m_e_MeV = 0.51099895
Z_MeV = m_e_MeV / (0.4 * np.pi)  # 0.40664005 MeV
Y2_target = 938.27208816 / Z_MeV  # m_p / Z ≈ 2307.3775
ell_P = np.sqrt(
    1.054571817e-34 * 6.67430e-11 / (299792458**3)
)  # 1.616255e-35 m
lambda_Z = (197.3269804e-15) / Z_MeV  # 4.85262e-13 m
Y3_target = ell_P / lambda_Z  # ell_P / lambda_Z ≈ 3.33069e-23

# 3. N.E.A. 理论基准计算
Y1_base = np.exp(5.0 * U_w) / R
Y2_base = f_g * (np.sqrt(R) / (1.0 + R)) * np.exp(0.5 * U_w)
Y3_base = ((1.0 + R) / f_g) * np.exp(-3.0 * U_w)

# 4. 对数残差 epsilon = ln(Y_base / Y_target)
eps1 = np.log(Y1_base / Y1_target)
eps2 = np.log(Y2_base / Y2_target)
eps3 = np.log(Y3_base / Y3_target)

print("=" * 72)
print("尺度骨架铁三角敏感度分析与对数全微分检验")
print("=" * 72)
print(f"理论基础值 vs 实验目标值:")
print(
    f"  Y1 (alpha_G^-1) : {Y1_base:.6e} vs {Y1_target:.6e} | 残差: {eps1*100:+.4f}%"
)
print(
    f"  Y2 (m_p / Z)    : {Y2_base:.4f} vs {Y2_target:.4f} | 残差: {eps2*100:+.4f}%"
)
print(
    f"  Y3 (ell_P/lam_Z): {Y3_base:.6e} vs {Y3_target:.6e} | 残差: {eps3*100:+.4f}%"
)
print("-" * 72)
print(f"指数恒等式闭环检验: eps1 + 2*eps2 + 2*eps3 = {eps1 + 2*eps2 + 2*eps3:.2e}")
print("结论: 3 个观测量残差严格满足代数闭环，仅有 2 个独立自由度。")
print("=" * 72)

# 5. 双通道修正结算
delta_p = -eps1 / U_w  # 修正指数 5
delta_g = -eps2  # 修正立体角因子 f_g

print(f"【双通道联合修正建议解】:")
print(f"  1. 主指数 5 的谱修正项 delta_p = {delta_p:+.8f} (有效指数 p = {5.0 + delta_p:.6f})")
print(f"  2. 立体角修正 delta_g = {delta_g:+.8f}")
print(f"     => 修正后有效立体角因子 f_g_new = {f_g * np.exp(delta_g):.6f}")
print(f"     => 对应二阶几何展开项: f_g = 1 + Delta/(4pi) - {abs(f_g - f_g * np.exp(delta_g)):.6f}")
print("=" * 72)