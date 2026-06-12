#!/usr/bin/env python3
"""
fit_tdep.py — Fit dipole t-dependence to Stepan's exponential.

Stepan:   dσ/dt ∝ exp(-b|t|),  b = 1.2 GeV⁻²
DIFFRAD:  dσ/dt ∝ 1/(mg² - t)^n    (t < 0, so mg²-t > 0)

Currently n=4, mg²=1.0.  Fit mg² and n to best match exp(-1.2|t|)
over the kinematic range |t| = 0 .. 5 GeV².
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

b_stepan = 1.2  # GeV^-2

# |t| grid  (t < 0 in physics, so |t| = -t)
mt = np.linspace(0.0, 5.0, 500)

# Stepan target (normalized to 1 at t=0)
f_stepan = np.exp(-b_stepan * mt)

# DIFFRAD dipole model:  f(|t|) = [mg² / (mg² + |t|)]^n
# (normalized to 1 at |t|=0)
def dipole(mt, mg2, n):
    return (mg2 / (mg2 + mt))**n

def chi2(params):
    mg2, n = params
    if mg2 < 0.05 or n < 0.5:
        return 1e12
    f = dipole(mt, mg2, n)
    # weight more at low |t| where most events are
    w = f_stepan  # weight by the target itself
    return np.sum(w * (f - f_stepan)**2)

# ── Scan mg² and n to visualize the landscape ────────────────────────────────
mg2_vals = np.linspace(0.3, 6.0, 100)
n_vals   = np.linspace(1.0, 8.0, 100)
chi2_map = np.zeros((len(n_vals), len(mg2_vals)))
for i, n in enumerate(n_vals):
    for j, mg2 in enumerate(mg2_vals):
        chi2_map[i, j] = chi2([mg2, n])

# ── Fit ──────────────────────────────────────────────────────────────────────
# Start from several initial points
best = None
for mg2_0 in [0.5, 1.0, 2.0, 3.0, 4.0]:
    for n_0 in [2.0, 3.0, 4.0, 5.0, 6.0]:
        res = minimize(chi2, [mg2_0, n_0], method='Nelder-Mead',
                       options=dict(xatol=1e-6, fatol=1e-10, maxiter=10000))
        if best is None or res.fun < best.fun:
            best = res

mg2_fit, n_fit = best.x
print(f"Best fit: mg² = {mg2_fit:.4f} GeV²,  n = {n_fit:.4f}")
print(f"  chi² = {best.fun:.6f}")
print(f"  b_eff(t=0) = n/mg² = {n_fit/mg2_fit:.3f} GeV⁻²")

# Also show current DIFFRAD and a few other combos
configs = [
    ("Current DIFFRAD", 1.0, 4.0),
    (f"Best fit", mg2_fit, n_fit),
    (f"n=4 best mg²", None, 4.0),   # fit mg² with n fixed at 4
    (f"n=3 best mg²", None, 3.0),   # fit mg² with n fixed at 3
    (f"n=2 best mg²", None, 2.0),   # fit mg² with n fixed at 2
]

# For fixed-n cases, find best mg²
for i, (label, mg2, n) in enumerate(configs):
    if mg2 is None:
        res_n = minimize(lambda x: chi2([x[0], n]), [2.0], method='Nelder-Mead',
                         options=dict(xatol=1e-6))
        configs[i] = (label, res_n.x[0], n)

print(f"\n{'Label':<20s} {'mg²':>8s} {'n':>6s} {'b_eff(0)':>10s} {'chi²':>12s}")
print("-" * 60)
for label, mg2, n in configs:
    c = chi2([mg2, n])
    print(f"{label:<20s} {mg2:8.4f} {n:6.2f} {n/mg2:10.3f} {c:12.6f}")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle(r'Fitting $1/(m_g^2 - t)^n$ to $e^{-1.2|t|}$  (Stepan)',
             fontsize=14, fontweight='bold')

# (a) Shape comparison
ax = axes[0]
ax.plot(mt, f_stepan, 'r-', lw=3, label=r'Stepan: $e^{-1.2|t|}$')
colors = ['gray', 'blue', 'green', 'orange', 'purple']
for (label, mg2, n), col in zip(configs, colors):
    ls = '--' if 'Current' in label else '-'
    lw = 1.5 if 'Current' in label else 2.0
    ax.plot(mt, dipole(mt, mg2, n), color=col, ls=ls, lw=lw,
            label=f'{label}: $m_g^2$={mg2:.2f}, n={n:.2f}')
ax.set_xlabel(r'$|t|$ (GeV$^2$)', fontsize=12)
ax.set_ylabel(r'd$\sigma$/d$t$ (normalized)', fontsize=12)
ax.set_title('Linear scale')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (b) Log scale
ax = axes[1]
ax.semilogy(mt, f_stepan, 'r-', lw=3, label=r'Stepan: $e^{-1.2|t|}$')
for (label, mg2, n), col in zip(configs, colors):
    ls = '--' if 'Current' in label else '-'
    lw = 1.5 if 'Current' in label else 2.0
    ax.semilogy(mt, dipole(mt, mg2, n), color=col, ls=ls, lw=lw,
                label=f'{label}')
ax.set_xlabel(r'$|t|$ (GeV$^2$)', fontsize=12)
ax.set_ylabel(r'd$\sigma$/d$t$ (normalized)', fontsize=12)
ax.set_title('Log scale')
ax.set_ylim(1e-4, 2)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (c) chi² map
ax = axes[2]
levels = np.logspace(np.log10(best.fun + 0.001), np.log10(chi2_map.max()), 30)
cs = ax.contourf(mg2_vals, n_vals, np.log10(chi2_map + 1e-10), levels=30, cmap='viridis')
plt.colorbar(cs, ax=ax, label=r'log$_{10}(\chi^2)$')
ax.plot(mg2_fit, n_fit, 'r*', ms=15, label=f'Best: ({mg2_fit:.2f}, {n_fit:.2f})')
ax.plot(1.0, 4.0, 'wo', ms=10, mew=2, label='Current (1.0, 4)')
for label, mg2, n in configs[2:]:
    ax.plot(mg2, n, 'w^', ms=8, mew=1.5)
ax.set_xlabel(r'$m_g^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel(r'$n$ (power)', fontsize=12)
ax.set_title(r'$\chi^2$ landscape')
ax.legend(fontsize=9)

fig.tight_layout()
fig.savefig('fit_tdep_result.png', dpi=150)
print(f"\nPlot saved: fit_tdep_result.png")

import os; os.system('open fit_tdep_result.png')
