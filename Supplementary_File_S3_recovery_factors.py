"""Supplementary File S3: recovery-factor, equivalent shut-in time and paired-difference
calculations reported in Section 4.5 and Table 7. Inputs are the pooled slopes (Table 5) and
per-well anchored gradients (Table 3) of the manuscript; no new BHT data are used."""
import numpy as np
from scipy.special import exp1
from scipy.optimize import brentq
rw, alpha, tc, dt_eff, tau, h = 0.108, 1.0e-6, 4.0, 24.3, 13.6, 3600.0
def R_horner(tc, dt):
    tc, dt = tc*h, dt*h
    return np.log((tc+dt)/dt)/exp1(rw**2/(4*alpha*tc))
def R_line(tc, dt):
    tc, dt = tc*h, dt*h
    return (exp1(rw**2/(4*alpha*(tc+dt))) - exp1(rw**2/(4*alpha*dt)))/exp1(rw**2/(4*alpha*tc))
G0, A0 = 25.35, 19.95   # uncorrected free-intercept slope and anchored gradient (degC/km)
free = {"110%":27.88, "Horner":29.00, "Line source":28.87, "Time-dependent":30.50}
anch = {"110%":22.73, "Horner":22.82, "Line source":22.72, "Time-dependent":24.00}
se   = {"Measured":(4.40,25.35),"110%":(4.84,27.88),"Horner":(5.03,29.0),"Line source":(5.01,28.87),"Time-dependent":(5.29,30.5)}
for k in free:
    a = free[k]/G0
    print(f"{k:15s} a={a:.4f} R={1-1/a:.4f}  anchored factor={anch[k]/A0:.4f}")
for k,(s,g) in se.items(): print(f"{k:15s} SE/slope={100*s/g:.2f}%  t={g/s:.2f}")
print("Model R (tc=4 h, dt=24.3 h): Horner %.4f, line source %.4f, exponential %.4f"
      % (R_horner(tc,dt_eff), R_line(tc,dt_eff), np.exp(-dt_eff/tau)))
a110 = 1.1; R110 = 1-1/a110
print("Horner dt for a=1.1: %.1f h; line-source dt: %.1f h" %
      (brentq(lambda d:R_horner(tc,d)-R110,1,500), brentq(lambda d:R_line(tc,d)-R110,1,500)))
print("Horner dt for anchored 110%% factor: %.1f h" %
      brentq(lambda d:R_horner(tc,d)-(1-A0/anch["110%"]),1,500))
print("Offset of Eq. 3 (0.1*Ts = 2.7 degC) in anchored gradient, z = 2.90..4.44 km:",
      [round(2.7/z,2) for z in (2.90,3.38,3.92,4.44)])
print("Eq. 3 applied to degF / K adds (degC):", round((1.1*32-32)/1.8,2), round(0.1*273.15,2))
print("Ts +/-2 degC changes the anchored gradient by", round(2/4.44,2), "-", round(2/2.90,2), "degC/km")
