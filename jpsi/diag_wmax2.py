#!/usr/bin/env python3
"""
Diagnostic #2: compare weight variance with uniform xB vs log(xB) sampling.
"""
import numpy as np

alpha = 1.0/137.036
barn  = 0.389379e6
Mp, Mp2 = 0.938272, 0.938272**2
me2   = 0.000511**2
Mjp, Mjp2 = 3.0969, 3.0969**2
Wth, Wth2 = Mp + Mjp, (Mp + Mjp)**2
pi = np.pi

alf1, alf2, alf3 = 400.0, 1.0, 0.32
nuT, mg2, cR = 3.0, 1.0, 0.4

q2min, q2max = 0.001, 0.25
xbmin, xbmax = 0.00001, 0.016
tmin, tmax = 0.0, 7.0
ymin, ymax = 0.79, 0.99
Wmin = 4.05
Ebeam = 10.6
s = 2*Mp*Ebeam + Mp2

def sigma_T_jpsi(q2, xB, t):
    if xB <= 0 or xB >= 1 or q2 <= 0: return 0.0
    w2 = Mp2 + q2*(1-xB)/xB
    if w2 <= Wth2: return 0.0
    cT = alf1 * (1 - Wth2/w2)**alf2 * np.sqrt(w2)**alf3
    sigT = cT / (1 + q2/Mjp2)**nuT
    return sigT * 3*mg2**3 / (mg2 - t)**4

def sigma_L_jpsi(q2, xB, t):
    return cR * (q2/Mjp2) * sigma_T_jpsi(q2, xB, t)

def bornin(q2, xB, y, t):
    gamma2 = 4*Mp2*xB**2/q2
    eps_num = 1 - y - 0.25*y**2*gamma2
    eps_den = 1 - y + 0.5*y**2 + 0.25*y**2*gamma2
    if eps_den <= 0 or eps_num <= 0: return 0.0
    eps = eps_num / eps_den
    sT = sigma_T_jpsi(q2, xB, t)
    sL = sigma_L_jpsi(q2, xB, t)
    return alpha*barn/(2*pi) * y**2/(1-eps) * (1-xB)/(xB*q2) * (sT + eps*sL)

def sample_uniform_xb(rng, N):
    r1, r2, r3 = rng.random(N), rng.random(N), rng.random(N)
    lq2min, lq2max = np.log(q2min), np.log(q2max)
    q2_arr = np.exp(lq2min + r1*(lq2max - lq2min))
    wjacq2 = q2_arr * (lq2max - lq2min)
    xb_arr = xbmin + r2*(xbmax - xbmin)
    wjacxb = xbmax - xbmin
    y_arr = q2_arr / (s * xb_arr)
    w2_arr = Mp2 + s*y_arr - q2_arr
    tdmin_u, tdmax_u = -tmax, -tmin
    absA = (mg2 - tdmin_u)**(-3)
    absB = (mg2 - tdmax_u)**(-3)
    abst = absB - absA
    ranexp = absA + r3*abst
    t_arr = mg2 - ranexp**(-1.0/3.0)
    wjact = abst * (mg2 - t_arr)**4 / 3.0
    sg_born = wjacq2 * wjacxb * wjact * 2*pi
    valid = ((xb_arr > 0) & (xb_arr < 1) &
             (y_arr > ymin) & (y_arr < ymax) & (y_arr < 1) &
             (w2_arr > Wmin**2) & (w2_arr > Wth2))
    sib = np.zeros(N)
    for i in np.where(valid)[0]:
        sib[i] = bornin(q2_arr[i], xb_arr[i], y_arr[i], t_arr[i])
    weight = sg_born * sib
    weight[~valid] = 0.0
    return weight, valid

def sample_log_xb(rng, N):
    r1, r2, r3 = rng.random(N), rng.random(N), rng.random(N)
    lq2min, lq2max = np.log(q2min), np.log(q2max)
    q2_arr = np.exp(lq2min + r1*(lq2max - lq2min))
    wjacq2 = q2_arr * (lq2max - lq2min)
    # Log sampling for xB
    lxbmin, lxbmax = np.log(xbmin), np.log(xbmax)
    xb_arr = np.exp(lxbmin + r2*(lxbmax - lxbmin))
    wjacxb = xb_arr * (lxbmax - lxbmin)
    y_arr = q2_arr / (s * xb_arr)
    w2_arr = Mp2 + s*y_arr - q2_arr
    tdmin_u, tdmax_u = -tmax, -tmin
    absA = (mg2 - tdmin_u)**(-3)
    absB = (mg2 - tdmax_u)**(-3)
    abst = absB - absA
    ranexp = absA + r3*abst
    t_arr = mg2 - ranexp**(-1.0/3.0)
    wjact = abst * (mg2 - t_arr)**4 / 3.0
    sg_born = wjacq2 * wjacxb * wjact * 2*pi
    valid = ((xb_arr > 0) & (xb_arr < 1) &
             (y_arr > ymin) & (y_arr < ymax) & (y_arr < 1) &
             (w2_arr > Wmin**2) & (w2_arr > Wth2))
    sib = np.zeros(N)
    for i in np.where(valid)[0]:
        sib[i] = bornin(q2_arr[i], xb_arr[i], y_arr[i], t_arr[i])
    weight = sg_born * sib
    weight[~valid] = 0.0
    return weight, valid

print("=== Comparing uniform xB vs log(xB) sampling ===\n")

N = 500000
print(f"--- Uniform xB sampling (N={N}) ---")
ntrial = 50
wmaxes_uni = []
for trial in range(ntrial):
    w, v = sample_uniform_xb(np.random.default_rng(42 + trial*7919), 10000)
    wmaxes_uni.append(w.max())
wmaxes_uni = np.array(wmaxes_uni)

w_uni, v_uni = sample_uniform_xb(np.random.default_rng(42), N)
m_uni = v_uni & (w_uni > 0)
print(f"Valid events: {m_uni.sum()}/{N} = {100*m_uni.sum()/N:.1f}%")
print(f"Weight: median={np.median(w_uni[m_uni]):.3e}, max={w_uni.max():.3e}, max/median={w_uni.max()/np.median(w_uni[m_uni]):.0f}x")
print(f"wmax over {ntrial} trials: median={np.median(wmaxes_uni):.3e}, max/min={wmaxes_uni.max()/wmaxes_uni.min():.1f}x")

print(f"\n--- Log xB sampling (N={N}) ---")
wmaxes_log = []
for trial in range(ntrial):
    w, v = sample_log_xb(np.random.default_rng(42 + trial*7919), 10000)
    wmaxes_log.append(w.max())
wmaxes_log = np.array(wmaxes_log)

w_log, v_log = sample_log_xb(np.random.default_rng(42), N)
m_log = v_log & (w_log > 0)
print(f"Valid events: {m_log.sum()}/{N} = {100*m_log.sum()/N:.1f}%")
print(f"Weight: median={np.median(w_log[m_log]):.3e}, max={w_log.max():.3e}, max/median={w_log.max()/np.median(w_log[m_log]):.0f}x")
print(f"wmax over {ntrial} trials: median={np.median(wmaxes_log):.3e}, max/min={wmaxes_log.max()/wmaxes_log.min():.1f}x")

# Effective efficiency comparison
wmax_true_uni = w_uni[m_uni].max()
wmax_true_log = w_log[m_log].max()
eff_uni = w_uni[m_uni].mean() / wmax_true_uni
eff_log = w_log[m_log].mean() / wmax_true_log
print(f"\n--- Effective efficiency (mean_weight / max_weight) ---")
print(f"Uniform xB: eff = {eff_uni:.4e}  (x{eff_log/eff_uni:.1f} worse)")
print(f"Log     xB: eff = {eff_log:.4e}")
print(f"\nImprovement factor: {eff_log/eff_uni:.1f}x faster generation with log(xB)")

# Also check: what about log sampling for BOTH Q2 and xB
# combined with the 1/(xB*Q2) factor in the flux?
print(f"\n--- Weight percentiles comparison ---")
print(f"{'percentile':>12} {'uniform xB':>14} {'log xB':>14} {'ratio':>8}")
for p in [50, 90, 95, 99, 99.9, 100]:
    v1 = np.percentile(w_uni[m_uni], p) if p < 100 else w_uni[m_uni].max()
    v2 = np.percentile(w_log[m_log], p) if p < 100 else w_log[m_log].max()
    print(f"{p:12.1f} {v1:14.3e} {v2:14.3e} {v1/v2:8.1f}")
