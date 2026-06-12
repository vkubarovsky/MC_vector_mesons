#!/usr/bin/env python3
"""
Diagnostic: map the weight landscape for J/psi Born generation
to understand why wmax varies 7x across jobs.

Replicates the Fortran sampling + cross-section in pure Python.
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
alpha = 1.0/137.036
barn  = 0.389379e6   # nb/GeV^2
Mp    = 0.938272
Mp2   = Mp**2
me2   = 0.000511**2
Mjp   = 3.0969
Mjp2  = Mjp**2
Wth   = Mp + Mjp
Wth2  = Wth**2
pi    = np.pi

# J/psi cross-section parameters
alf1, alf2, alf3 = 400.0, 1.0, 0.32
nuT  = 3.0
mg2  = 1.0
cR   = 0.4

# Input ranges
q2min, q2max = 0.001, 0.25
xbmin, xbmax = 0.00001, 0.016
tmin,  tmax  = 0.0, 7.0
ymin,  ymax  = 0.79, 0.99
Wmin  = 4.05
Ebeam = 10.6
s     = 2*Mp*Ebeam + Mp2

def sigma_T_jpsi(q2, xB, t):
    if xB <= 0 or xB >= 1 or q2 <= 0:
        return 0.0
    w2 = Mp2 + q2*(1-xB)/xB
    if w2 <= Wth2:
        return 0.0
    cT = alf1 * (1 - Wth2/w2)**alf2 * np.sqrt(w2)**alf3
    sigT = cT / (1 + q2/Mjp2)**nuT
    return sigT * 3*mg2**3 / (mg2 - t)**4

def sigma_L_jpsi(q2, xB, t):
    sT = sigma_T_jpsi(q2, xB, t)
    return cR * (q2/Mjp2) * sT

def bornin(q2, xB, y, t):
    """d3sigma/(dxB dQ2 dt) in nb/GeV^4"""
    gamma2 = 4*Mp2*xB**2/q2
    eps_num = 1 - y - 0.25*y**2*gamma2
    eps_den = 1 - y + 0.5*y**2 + 0.25*y**2*gamma2
    if eps_den <= 0 or eps_num <= 0:
        return 0.0
    eps = eps_num / eps_den
    sT = sigma_T_jpsi(q2, xB, t)
    sL = sigma_L_jpsi(q2, xB, t)
    return alpha*barn/(2*pi) * y**2/(1-eps) * (1-xB)/(xB*q2) * (sT + eps*sL)

def sample_and_weight(rng, N):
    """Replicate the Fortran sampling, return (Q2, xB, t, weight) arrays."""
    r1 = rng.random(N)
    r2 = rng.random(N)
    r3 = rng.random(N)

    # Q2: log-sampling
    lq2min, lq2max = np.log(q2min), np.log(q2max)
    q2_arr = np.exp(lq2min + r1*(lq2max - lq2min))
    wjacq2 = q2_arr * (lq2max - lq2min)

    # xB: uniform
    xb_arr = xbmin + r2*(xbmax - xbmin)
    wjacxb = xbmax - xbmin

    # Derived: y = Q2/(s*xB)
    y_arr = q2_arr / (s * xb_arr)

    # W2
    w2_arr = Mp2 + s*y_arr - q2_arr

    # t: dipole sampling  1/(mg2-t)^4
    tdmin_u = -tmax   # = -7
    tdmax_u = -tmin   # = 0
    absA = (mg2 - tdmin_u)**(-3)   # (1+7)^{-3}
    absB = (mg2 - tdmax_u)**(-3)   # (1-0)^{-3} = 1
    abst = absB - absA
    ranexp = absA + r3*abst
    t_arr = mg2 - ranexp**(-1.0/3.0)
    wjact = abst * (mg2 - t_arr)**4 / 3.0

    sg_born = wjacq2 * wjacxb * wjact * 2*pi

    # Apply kinematic cuts (matching Fortran)
    valid = ((xb_arr > 0) & (xb_arr < 1) &
             (y_arr > ymin) & (y_arr < ymax) & (y_arr < 1) &
             (w2_arr > Wmin**2) & (w2_arr > Wth2))

    # Compute cross-section for valid events
    sib = np.zeros(N)
    for i in np.where(valid)[0]:
        sib[i] = bornin(q2_arr[i], xb_arr[i], y_arr[i], t_arr[i])

    weight = sg_born * sib
    weight[~valid] = 0.0

    return q2_arr, xb_arr, y_arr, t_arr, w2_arr, weight, valid

# Run large sample
print("Sampling 2,000,000 points...")
rng = np.random.default_rng(42)
q2, xb, y, t, w2, weight, valid = sample_and_weight(rng, 2_000_000)

# Only look at valid events
mask = valid & (weight > 0)
q2v  = q2[mask]
xbv  = xb[mask]
yv   = y[mask]
tv   = t[mask]
w2v  = w2[mask]
wv   = weight[mask]

print(f"Valid events: {mask.sum()} / {len(weight)}")
print(f"Weight range: {wv.min():.3e} to {wv.max():.3e}")
print(f"Weight mean:  {wv.mean():.3e}")
print(f"Weight median:{np.median(wv):.3e}")

# Look at the distribution of max-of-10000 samples
print("\n--- Simulating 100 jobs of 10,000 warmup points ---")
wmaxes = []
for trial in range(100):
    rng_trial = np.random.default_rng(333522 + trial*100003)
    _, _, _, _, _, w_trial, v_trial = sample_and_weight(rng_trial, 10000)
    wmaxes.append(w_trial.max())
wmaxes = np.array(wmaxes)
print(f"wmax range: {wmaxes.min():.3e} to {wmaxes.max():.3e}")
print(f"wmax mean:  {wmaxes.mean():.3e}")
print(f"wmax median:{np.median(wmaxes):.3e}")
print(f"wmax std:   {wmaxes.std():.3e}")
print(f"wmax ratio max/median: {wmaxes.max()/np.median(wmaxes):.2f}")

# Find the top 20 weights and their kinematics
top_idx = np.argsort(wv)[-20:][::-1]
print(f"\n--- Top 20 weights ---")
print(f"{'rank':>4} {'weight':>12} {'Q2':>10} {'xB':>12} {'y':>8} {'t':>8} {'W':>8}")
for rank, idx in enumerate(top_idx):
    print(f"{rank+1:4d} {wv[idx]:12.3e} {q2v[idx]:10.6f} {xbv[idx]:12.8f} "
          f"{yv[idx]:8.4f} {tv[idx]:8.4f} {np.sqrt(w2v[idx]):8.4f}")

# Percentile analysis
print("\n--- Weight percentiles ---")
for p in [50, 90, 95, 99, 99.5, 99.9, 99.95, 99.99, 100]:
    val = np.percentile(wv, p) if p < 100 else wv.max()
    print(f"  {p:7.2f}%: {val:.3e}")

# ---- Plots ----
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# 1. Weight distribution (log scale)
ax = axes[0, 0]
logw = np.log10(wv)
ax.hist(logw, bins=100, color='steelblue', edgecolor='none')
ax.set_xlabel('log10(weight)')
ax.set_ylabel('counts')
ax.set_title('Weight distribution')
ax.axvline(np.log10(np.median(wmaxes)), color='red', ls='--', label=f'median wmax={np.median(wmaxes):.2e}')
ax.axvline(np.log10(wmaxes.max()), color='red', ls='-', label=f'max wmax={wmaxes.max():.2e}')
ax.legend(fontsize=8)

# 2. Weight vs Q2
ax = axes[0, 1]
ax.scatter(q2v, wv, s=0.1, alpha=0.1, c='steelblue')
ax.set_xlabel('Q2 (GeV^2)')
ax.set_ylabel('weight')
ax.set_yscale('log')
ax.set_title('Weight vs Q2')

# 3. Weight vs xB
ax = axes[0, 2]
ax.scatter(xbv, wv, s=0.1, alpha=0.1, c='steelblue')
ax.set_xlabel('xB')
ax.set_ylabel('weight')
ax.set_yscale('log')
ax.set_title('Weight vs xB')

# 4. Weight vs t
ax = axes[1, 0]
ax.scatter(tv, wv, s=0.1, alpha=0.1, c='steelblue')
ax.set_xlabel('t (GeV^2)')
ax.set_ylabel('weight')
ax.set_yscale('log')
ax.set_title('Weight vs t')

# 5. Weight vs y
ax = axes[1, 1]
ax.scatter(yv, wv, s=0.1, alpha=0.1, c='steelblue')
ax.set_xlabel('y')
ax.set_ylabel('weight')
ax.set_yscale('log')
ax.set_title('Weight vs y')

# 6. Weight vs W
ax = axes[1, 2]
ax.scatter(np.sqrt(w2v), wv, s=0.1, alpha=0.1, c='steelblue')
ax.set_xlabel('W (GeV)')
ax.set_ylabel('weight')
ax.set_yscale('log')
ax.set_title('Weight vs W')

fig.suptitle('J/psi weight landscape diagnostic', fontsize=14)
fig.tight_layout()
fig.savefig('diag_wmax.png', dpi=150)
print("\nSaved: diag_wmax.png")

# ---- Distribution of wmax across warmup sizes ----
fig2, ax2 = plt.subplots(figsize=(8, 5))
for nwarm_test in [1000, 5000, 10000, 50000, 100000]:
    wmaxes_test = []
    for trial in range(50):
        rng_t = np.random.default_rng(42 + trial*7919 + nwarm_test)
        _, _, _, _, _, w_t, _ = sample_and_weight(rng_t, nwarm_test)
        wmaxes_test.append(w_t.max())
    wmaxes_test = np.array(wmaxes_test)
    ratio = wmaxes_test.max() / wmaxes_test.min()
    ax2.hist(wmaxes_test/1e8, bins=20, alpha=0.5, label=f'nwarm={nwarm_test} (max/min={ratio:.1f}x)')

ax2.set_xlabel('wmax (x 10^8)')
ax2.set_ylabel('count (50 trials)')
ax2.set_title('wmax distribution vs warmup size')
ax2.legend()
fig2.tight_layout()
fig2.savefig('diag_wmax_vs_nwarm.png', dpi=150)
print("Saved: diag_wmax_vs_nwarm.png")
