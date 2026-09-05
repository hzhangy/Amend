import numpy as np

print("="*70)
print("  N.E.A. 样本外检验：电弱扇区 sin^2(theta_W) 与 32pi^2 因子")
print("="*70)

# 1. N.E.A. 拓扑常数
pi = np.pi
sqrt3 = np.sqrt(3.0)
R = 1.0 / (1.0 + pi)
Delta = 1.0 - sqrt3 / 2.0
delta = Delta / (4.0 * pi)

# 2. sin^2(theta_W) 的 N.E.A. 原始公式
# 理论值 = R - delta
s_theo = R - delta

# 3. 实验值 (PDG 2022/2023, MS-bar scheme at M_Z)
# 有效温伯格角 sin^2(theta_eff) ≈ 0.23122
s_exp = 0.23122 

# 4. 计算绝对残差
residual = s_exp - s_theo

# 5. 测试候选修正因子：Delta / (32 * pi^2)
candidate_factor = Delta / (32.0 * pi**2)

# 6. 比例计算
ratio = residual / candidate_factor

print(f"[N.E.A. 理论值] sin^2(theta_W) = R - delta      = {s_theo:.6f}")
print(f"[PDG 实验值]    sin^2(theta_W) (MS-bar)       = {s_exp:.6f}")
print("-" * 70)
print(f"[绝对残差]      Exp - Theo                    = {residual:.6f}  ({residual*1e6:.1f} ppm)")
print(f"[候选因子]      Delta / (32 * pi^2)           = {candidate_factor:.6f}  ({candidate_factor*1e6:.1f} ppm)")
print("="*70)
print(f"\n【核心检验：残差 / 候选因子】")
print(f"  Ratio = {ratio:.4f}")
print(f"  (如果 Ratio 接近 1.0，则证明 32pi^2 是跨扇区普适因子！)")
print("="*70)