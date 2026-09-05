import numpy as np

print("="*75)
print("  N.E.A. 铁三角高阶修正分配验证 (策略1: 修正归入空间几何 ell_P)")
print("="*75)

# 1. 实验基准值 (基于 CODATA 2018 与 N.E.A. 标度 Z)
# 为保持与上一轮审计一致，直接使用已确认的原始偏差 ppm 值
err_alpha_raw_ppm = -1203.0  # alpha_G^-1 原始偏小 1203 ppm
err_ell_raw_ppm   =  430.0   # ell_P/lambda_Z 原始偏大 430 ppm
err_m_raw_ppm     =  170.0   # m_p/Z 原始偏大 170 ppm

# 2. N.E.A. 拓扑修正因子
Delta = 1.0 - np.sqrt(3.0) / 2.0
C = 1.0 + (Delta**3) / 2.0   # C ≈ 1.00120237

# 修正量 (ppm)
C_ppm = (C - 1.0) * 1e6      # +1202.37 ppm
C_inv_half_ppm = (C**(-0.5) - 1.0) * 1e6  # -601.18 ppm

# 3. 策略 1 分配计算
# alpha_G^-1 吸收全部 C (即 +1202.37 ppm)
err_alpha_corr = err_alpha_raw_ppm + C_ppm

# ell_P/lambda_Z 吸收 C^(-1/2) (即 -601.18 ppm)
err_ell_corr = err_ell_raw_ppm + C_inv_half_ppm

# m_p/Z 保持不变
err_m_corr = err_m_raw_ppm

# 4. 打印审计报告
print(f"{'参数':<18} | {'原始偏差 (ppm)':>15} | {'修正后偏差 (ppm)':>18}")
print("-" * 75)
print(f"{'alpha_G^-1':<18} | {err_alpha_raw_ppm:>15.1f} | {err_alpha_corr:>18.1f}")
print(f"{'ell_P / lambda_Z':<18} | {err_ell_raw_ppm:>15.1f} | {err_ell_corr:>18.1f}")
print(f"{'m_p / Z':<18} | {err_m_raw_ppm:>15.1f} | {err_m_corr:>18.1f}")
print("="*75)

print("\n【闭合校验】:")
print(f"  修正后闭合残差 = -2 * (ell_P 偏差) - 2 * (m_p 偏差) - (alpha_G 偏差)")
closure_residual = -2 * err_ell_corr - 2 * err_m_corr - err_alpha_corr
print(f"  计算结果: {closure_residual:.2f} ppm (必须接近 0 以维持物理恒等式)")