import numpy as np
Delta = 1.0 - np.sqrt(3.0)/2.0
a_nea = 25*np.sqrt(3)*np.pi + 1.0
a_exp = 137.035999084
C3 = 1.0 + Delta**3/(32*np.pi**2)
a_corr = a_nea * C3
print(f"raw ppm      = {(a_nea-a_exp)/a_exp*1e6:+.3f}")
print(f"C3-1  (ppm)  = {(C3-1)*1e6:+.3f}")
print(f"corrected ppm= {(a_corr-a_exp)/a_exp*1e6:+.3f}")