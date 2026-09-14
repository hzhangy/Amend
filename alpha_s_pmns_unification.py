#!/usr/bin/env python3
"""
alpha_s_pmns_unification.py

从 α_s 的 1/28 出发，探索统一模式。

关键观察：
  - α_s 修正 ~ 1/28 = 1/E(K8) = 1/C(8,2)
  - m_μ 的 N = 36 = E(K9) = C(9,2)
  - m_0 的 9 = V(C8)+1

所有这些都涉及 C8 顶点数 8 的组合。

探索：
  1. α_s 修正是否也用"完全图边数"模式
  2. PMNS θ_13 是否有类似的几何表达
  3. 是否有统一公式
"""

import numpy as np
from itertools import combinations

pi = np.pi
alpha_inv = 137.035999084
Delta = 1 - np.sqrt(3)/2
R = 1/(1 + pi)

# ============================================================
# α_s
# ============================================================
alpha_s_pred = 2/(10*np.sqrt(3)) * (1 + R**2)
alpha_s_obs = 0.1179

print("=" * 90)
print("  α_s 的 1/28 探索")
print("=" * 90)
print()

print(f"  α_s 主项 = {alpha_s_pred:.6f}")
print(f"  α_s 观测 = {alpha_s_obs}")
print(f"  主项偏差 = {(alpha_s_pred - alpha_s_obs)/alpha_s_obs:+.4%}")
print()

# 需要的修正
delta_needed = 1 - alpha_s_obs / alpha_s_pred
print(f"  需要的修正 = 1 - {delta_needed:.6f} = 1/{1/delta_needed:.4f}")
print()

# 完全图边数
EKn = {n: n*(n-1)//2 for n in range(3, 20)}
print("  完全图边数 E(K_n)：")
for n in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
    print(f"    E(K_{n}) = {EKn[n]}")
print()

# 精确需要的分母
N_needed = 1 / delta_needed
print(f"  精确需要的分母 N = {N_needed:.4f}")
print()

# 匹配
print("  完全图匹配：")
for n, val in EKn.items():
    if abs(val - N_needed) / N_needed < 0.05:
        # 验证
        alpha_s_new = alpha_s_pred * (1 - 1/val)
        err = abs(alpha_s_new - alpha_s_obs) / alpha_s_obs
        print(f"    N = E(K_{n}) = {val}: α_s = {alpha_s_new:.6f}  (偏差 {err:+.4%})")
print()

# ============================================================
# α_s 修正的物理意义
# ============================================================
print("=" * 90)
print("  α_s 修正的物理意义")
print("=" * 90)
print()

print("  E(K_8) = 28 = C(8,2) = 8×7/2")
print("  V(C8) = 8")
print()
print("  物理解释：")
print("    α_s 是强耦合常数")
print("    修正 = 1/C(V(C8), 2) = 1/C(8,2) = 1/28")
print("    = 'C8 顶点之间所有可能连接数的倒数'")
print()
print("  与 m_μ 的对比：")
print("    m_μ 修正 N = 36 = C(9,2) = C(V(C8)+1, 2)")
print("    α_s 修正 N = 28 = C(8,2) = C(V(C8), 2)")
print()
print("  统一形式：")
print("    修正率 = 1/C(n,2)")
print("    其中 n 是某个几何顶点数")
print()

# ============================================================
# PMNS θ_13
# ============================================================
print("=" * 90)
print("  PMNS θ_13 探索")
print("=" * 90)
print()

theta_13_obs = 8.58  # degrees
sin_theta_13_obs = np.sin(np.radians(theta_13_obs))

print(f"  θ_13 观测 = {theta_13_obs}°")
print(f"  sin(θ_13) = {sin_theta_13_obs:.6f}")
print()

# 系统搜索：用 Δ 的幂和 1/n 形式
print("  系统搜索（Δ 的幂 × 有理数）：")
print()

# 尝试 sin(θ_13) = Δ^p × (a/b)
candidates = []
for p in [0.5, 1.0, 1.1, 1.2, 1.5, 2.0]:
    for a in range(1, 11):
        for b in range(1, 21):
            val = (Delta ** p) * (a / b)
            err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
            if err < 0.005:
                candidates.append((f"Δ^{p} × {a}/{b}", val, err))

# 尝试 sin(θ_13) = R^p × (a/b)
for p in [0.5, 1.0, 1.1, 1.2, 1.5, 2.0]:
    for a in range(1, 11):
        for b in range(1, 21):
            val = (R ** p) * (a / b)
            err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
            if err < 0.005:
                candidates.append((f"R^{p} × {a}/{b}", val, err))

# 尝试 sin(θ_13) = Δ^p / (a/b) 的反向
for p in [0.5, 1.0, 1.1, 1.2, 1.5]:
    for a in range(1, 11):
        for b in range(1, 21):
            val = (Delta ** p) / (a / b)
            err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
            if err < 0.005:
                candidates.append((f"Δ^{p} / ({a}/{b})", val, err))

# 尝试 (a/b) × Δ + (c/d) × R
for a in range(1, 6):
    for b in range(1, 11):
        for c in range(1, 6):
            for d in range(1, 11):
                val = (a/b) * Delta + (c/d) * R
                err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
                if err < 0.002:
                    candidates.append((f"{a}/{b}×Δ + {c}/{d}×R", val, err))

# 排序并显示前 20 个
candidates.sort(key=lambda x: x[2])
print(f"  找到 {len(candidates)} 个候选（误差 < 0.5%）：")
print()
print(f"  {'表达式':<30} | {'值':>10} | {'偏差':>10}")
print("  " + "-" * 56)
for expr, val, err in candidates[:20]:
    print(f"  {expr:<30} | {val:>10.6f} | {err:>9.4%}")
print()

# ============================================================
# 用 1/n 形式表达 θ_13
# ============================================================
print("=" * 90)
print("  θ_13 的 1/n 形式")
print("=" * 90)
print()

# sin(θ_13) 可能 = 1/n + small correction
# 1/sin(θ_13) = 1/0.149 = 6.70
print(f"  1/sin(θ_13) = {1/sin_theta_13_obs:.4f}")
print()

# 候选 1/n
for n in range(3, 30):
    val = 1/n
    err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
    if err < 0.1:
        print(f"    1/{n} = {val:.4f}  (偏差 {err:+.2%})")
print()

# sin(θ_13) = 1/(2π + x)?
print(f"  1/sin(θ_13) = {1/sin_theta_13_obs:.4f}")
print(f"  2π = {2*pi:.4f}")
print(f"  1/sin(θ_13) - 2π = {1/sin_theta_13_obs - 2*pi:.4f}")
print()

# 1/sin(θ_13) 的可能几何形式
print(f"  1/sin(θ_13) 的几何候选：")
print(f"    π + gen = {pi + 3:.4f}")
print(f"    π + √3 = {pi + np.sqrt(3):.4f}")
print(f"    2π + Δ = {2*pi + Delta:.4f}")
print(f"    2π + 1/2 = {2*pi + 0.5:.4f}")
print(f"    V_octa + 1 = 7")  # 6+1
print(f"    V_C8 - 1 = 7")
print(f"    E(K4) + 1 = 7")
print()

# ============================================================
# 从 sin(θ_13) 探索
# ============================================================
print("=" * 90)
print("  sin(θ_13) 的几何形式")
print("=" * 90)
print()

print(f"  sin(θ_13) = {sin_theta_13_obs:.6f}")
print()

# 尝试：sin(θ_13) = √(Δ/n)
print(f"  √(Δ/n) 形式：")
for n in [4, 5, 6, 7, 8, 9, 10]:
    val = np.sqrt(Delta/n)
    err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
    marker = " ✓" if err < 0.05 else ""
    print(f"    √(Δ/{n}) = {val:.6f}  (偏差 {err:+.2%}){marker}")
print()

# 尝试：sin(θ_13) = Δ^(3/2)
print(f"  Δ 的幂：")
for p in [0.5, 1.0, 1.1, 1.2, 1.3, 1.5]:
    val = Delta ** p
    err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
    marker = " ✓" if err < 0.05 else ""
    print(f"    Δ^{p} = {val:.6f}  (偏差 {err:+.2%}){marker}")
print()

# ============================================================
# 尝试：结合 D 和 ± 的模式
# ============================================================
print("=" * 90)
print("  尝试：PMNS 角是否用 D 的模式")
print("=" * 90)
print()

# θ_12 已知 = arctan((1-Δ/2)/√2)
# θ_23 = 45° 
# θ_13 = ?

# 尝试 θ_13 用 1/N 形式
# sin(θ_13) = ε/N?
for N in [3, 5, 6, 8, 9, 27, 28, 36, 60]:
    val = (1/10) / N
    err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
    if err < 0.5:
        print(f"    ε/{N} = {val:.6f}  (偏差 {err:+.2%})")
print()

# 尝试 sin(θ_13) = Δ/N
print(f"  Δ/N 形式：")
for N in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    val = Delta/N
    err = abs(val - sin_theta_13_obs) / sin_theta_13_obs
    if err < 0.2:
        print(f"    Δ/{N} = {val:.6f}  (偏差 {err:+.2%})")
print()

# ============================================================
# 综合统一尝试
# ============================================================
print("=" * 90)
print("  综合统一尝试")
print("=" * 90)
print()

print("  观察：α_s 和 m_μ 的修正都涉及完全图边数")
print()

print("  统一形式假设：")
print("    修正 = 1/C(n,2)")
print()
print("  具体：")
print("    α_s 修正 = 1/E(K_8) = 1/28  →  n = 8 = V(C8)")
print("    m_μ 修正 N = E(K_9) = 36    →  n = 9 = V(C8) + 1")
print()

print("  问题：为什么 α_s 用 n=8，m_μ 用 n=9？")
print("    可能：α_s 是规范耦合（直接与 C8 空间），用 V(C8)")
print("         m_μ 是物质质量（涉及弱+代），用 V(C8)+1")
print()

# ============================================================
# 精确验证
# ============================================================
print("=" * 90)
print("  精确验证：α_s 用 E(K8)")
print("=" * 90)
print()

alpha_s_EK8 = alpha_s_pred * (1 - 1/EKn[8])
print(f"  α_s 主项：{alpha_s_pred:.6f}")
print(f"  修正：1 - 1/28 = {1 - 1/28:.6f}")
print(f"  α_s (修正)：{alpha_s_EK8:.6f}")
print(f"  α_s 观测：{alpha_s_obs}")
print(f"  偏差：{(alpha_s_EK8 - alpha_s_obs)/alpha_s_obs:+.4%}")
print()

# 改进：加二级修正？
print("  是否需要二级修正？")
target_ratio = alpha_s_obs / alpha_s_pred
print(f"  目标 ratio：{target_ratio:.6f}")
print(f"  一级修正：{1 - 1/28:.6f}")
print(f"  差：{target_ratio - (1 - 1/28):.6f}")
print()

# 尝试二级
if abs(target_ratio - (1 - 1/28)) > 1e-5:
    diff = 1 - 1/28 - target_ratio
    print(f"  需要额外修正：{diff:.6f}")
    print(f"  = 1/{1/diff:.2f}")
    print()

# ============================================================
# 总结
# ============================================================
print("=" * 90)
print("  总结")
print("=" * 90)
print()

print("  1. α_s 的 1/28 = 1/E(K8) = 1/C(V(C8), 2)")
print(f"     修正后偏差：{(alpha_s_EK8 - alpha_s_obs)/alpha_s_obs:+.4%}")
print()
print("  2. PMNS θ_13 的候选：")
if candidates:
    print(f"     最佳：{candidates[0][0]} = {candidates[0][1]:.6f}")
    print(f"     偏差：{candidates[0][2]:.4%}")
print()
print("  3. 统一模式：")
print("     修正率 = 1/C(n, 2)")
print("     α_s: n = V(C8) = 8")
print("     m_μ: n = V(C8) + 1 = 9")
print()

print("  4. 仍待攻：")
print("     - 为什么 α_s 用 n=8，m_μ 用 n=9")
print("     - PMNS θ_13 的第一性来源")
print()