#!/usr/bin/env python3
"""
fit_Q2dep_M2.py — Fit M² in the Q² dependence with νT=3 fixed.

DIFFRAD:  σ_T(Q²) ∝ 1/(1 + Q²/M²)^νT
Current:  νT = 3,  M² = m²_J/ψ = 9.59 GeV²

Now let M² float (keep νT=3).

Usage:
    python fit_Q2dep_M2.py -diffrad born_events.lund -mariana Mariana.lund
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar
import argparse, os
from lund_reader import read_lund_events

parser = argparse.ArgumentParser()
parser.add_argument('-diffrad', required=True)
parser.add_argument('-mariana', required=True)
parser.add_argument('-nev', type=int, default=0)
args = parser.parse_args()

mJsq = 3.0969**2  # 9.59 GeV²

def load_Q2(fname, nev):
    print(f"Reading {fname} ...")
    evts = read_lund_events(fname, max_ev=nev)
    Q2 = np.array([e['Q2'] for e in evts])
    print(f"  {len(evts)} events,  Q² range: [{Q2.min():.4f}, {Q2.max():.4f}]")
    return Q2

dQ2 = load_Q2(args.diffrad, args.nev)
mQ2 = load_Q2(args.mariana, args.nev)

nuT_fix = 3.0
M2_0    = mJsq  # current: 9.59

def sigma_Q2(Q2, M2):
    return 1.0 / (1.0 + Q2/M2)**nuT_fix

sig0 = sigma_Q2(dQ2, M2_0)

bins_Q2 = np.linspace(0.001, 0.25, 35)
h_mar, _ = np.histogram(mQ2, bins=bins_Q2)
h_mar = h_mar.astype(float)
norm_mar = h_mar.sum()
if norm_mar > 0:
    h_mar /= norm_mar

def chi2_Q2(M2):
    if M2 < 0.001:
        return 1e12
    sig_new = sigma_Q2(dQ2, M2)
    w = sig_new / sig0
    h_dif, _ = np.histogram(dQ2, bins=bins_Q2, weights=w)
    s = h_dif.sum()
    if s > 0:
        h_dif /= s
    ok = h_mar > 0
    diff = h_dif[ok] - h_mar[ok]
    err2 = h_mar[ok] / norm_mar + 1e-6
    return np.sum(diff**2 / err2)

# ── Scan M² ──────────────────────────────────────────────────────────────────
M2_scan = np.logspace(-2, 2, 300)  # 0.01 to 100
chi2_scan = np.array([chi2_Q2(m) for m in M2_scan])

# ── Minimize ─────────────────────────────────────────────────────────────────
res = minimize_scalar(chi2_Q2, bounds=(0.001, 100.0), method='bounded')
M2_fit = res.x

print(f"\nFixed νT = {nuT_fix}")
print(f"Best fit:  M² = {M2_fit:.4f} GeV²  (M = {np.sqrt(M2_fit):.4f} GeV)")
print(f"  χ²(Q²) = {res.fun:.2f}")
print(f"  Current M²=mJ²={M2_0:.2f}:  χ²(Q²) = {chi2_Q2(M2_0):.2f}")

# Also try fitting both M² and νT
def chi2_both(params):
    M2, nuT = params
    if M2 < 0.001 or nuT < 0.5:
        return 1e12
    sig_new = 1.0 / (1.0 + dQ2/M2)**nuT
    w = sig_new / sig0
    h_dif, _ = np.histogram(dQ2, bins=bins_Q2, weights=w)
    s = h_dif.sum()
    if s > 0:
        h_dif /= s
    ok = h_mar > 0
    diff = h_dif[ok] - h_mar[ok]
    err2 = h_mar[ok] / norm_mar + 1e-6
    return np.sum(diff**2 / err2)

best2 = None
for m2_0 in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
    for nuT_0 in [1.0, 2.0, 3.0, 4.0, 5.0]:
        r = minimize(chi2_both, [m2_0, nuT_0], method='Nelder-Mead',
                     options=dict(xatol=1e-6, fatol=1e-6, maxiter=20000))
        if best2 is None or r.fun < best2.fun:
            best2 = r

M2_fit2, nuT_fit2 = best2.x
print(f"\nFree fit (M² and νT):")
print(f"  M² = {M2_fit2:.4f} GeV²  (M = {np.sqrt(max(M2_fit2,0)):.4f} GeV)")
print(f"  νT = {nuT_fit2:.2f}")
print(f"  χ²(Q²) = {best2.fun:.2f}")

# Table
print(f"\n{'M²':>8s} {'νT':>6s} {'M (GeV)':>10s} {'χ²(Q²)':>10s}")
print("-" * 40)
for M2_try in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, mJsq]:
    c = chi2_Q2(M2_try)
    tag = "  ← current" if abs(M2_try - mJsq) < 0.01 else ""
    print(f"{M2_try:8.3f} {nuT_fix:6.1f} {np.sqrt(M2_try):10.4f} {c:10.2f}{tag}")
print(f"{M2_fit:8.4f} {nuT_fix:6.1f} {np.sqrt(M2_fit):10.4f} {res.fun:10.2f}  ← best (νT=3)")
print(f"{M2_fit2:8.4f} {nuT_fit2:6.2f} {np.sqrt(max(M2_fit2,0)):10.4f} {best2.fun:10.2f}  ← best (free)")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle(r'Fitting $Q^2$: $1/(1+Q^2/M^2)^{\nu_T}$ with $\nu_T=3$ fixed, $M^2$ free',
             fontsize=14, fontweight='bold')

# (a) Q² distributions linear
ax = axes[0]
ax.hist(mQ2, bins=bins_Q2, density=True, histtype='stepfilled',
        color='lightsalmon', edgecolor='darkred', alpha=0.7, label='Mariana')
ax.hist(dQ2, bins=bins_Q2, density=True, histtype='step',
        color='navy', lw=1.5, ls='--', label=f'DIFFRAD orig (M²={M2_0:.1f})')

w_fit = sigma_Q2(dQ2, M2_fit) / sig0
ax.hist(dQ2, bins=bins_Q2, weights=w_fit, density=True, histtype='step',
        color='blue', lw=2.5, label=f'νT=3, M²={M2_fit:.3f}')

# Free fit
w_fit2 = (1.0/(1.0 + dQ2/M2_fit2)**nuT_fit2) / sig0
ax.hist(dQ2, bins=bins_Q2, weights=w_fit2, density=True, histtype='step',
        color='red', lw=2, ls='-.', label=f'νT={nuT_fit2:.1f}, M²={M2_fit2:.3f}')

for M2_show, col in [(0.05, 'green'), (0.2, 'orange'), (1.0, 'purple')]:
    w_s = sigma_Q2(dQ2, M2_show) / sig0
    ax.hist(dQ2, bins=bins_Q2, weights=w_s, density=True, histtype='step',
            color=col, lw=1.5, label=f'M²={M2_show}')

ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel('Normalized', fontsize=10)
ax.set_title(r'$Q^2$ distribution')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# (b) log scale
ax = axes[1]
ax.hist(mQ2, bins=bins_Q2, density=True, histtype='stepfilled',
        color='lightsalmon', edgecolor='darkred', alpha=0.7, label='Mariana')
ax.hist(dQ2, bins=bins_Q2, density=True, histtype='step',
        color='navy', lw=1.5, ls='--', label=f'orig (M²={M2_0:.1f})')
ax.hist(dQ2, bins=bins_Q2, weights=w_fit, density=True, histtype='step',
        color='blue', lw=2.5, label=f'νT=3, M²={M2_fit:.3f}')
ax.hist(dQ2, bins=bins_Q2, weights=w_fit2, density=True, histtype='step',
        color='red', lw=2, ls='-.', label=f'νT={nuT_fit2:.1f}, M²={M2_fit2:.3f}')
ax.set_yscale('log')
ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel('Normalized', fontsize=10)
ax.set_title(r'$Q^2$ distribution (log)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (c) χ² vs M²
ax = axes[2]
ax.semilogx(M2_scan, chi2_scan, 'b-', lw=2)
ax.axvline(M2_fit, color='red', ls='--', lw=1.5, label=f'Best: M²={M2_fit:.3f}')
ax.axvline(M2_0, color='gray', ls=':', lw=1.5, label=f'Current: M²={M2_0:.1f}')
ax.set_xlabel(r'$M^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel(r'$\chi^2(Q^2)$', fontsize=12)
ax.set_title(r'$\chi^2$ vs $M^2$ ($\nu_T=3$ fixed)')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('fit_Q2dep_M2_result.png', dpi=150)
print(f"\nPlot saved: fit_Q2dep_M2_result.png")
os.system('open fit_Q2dep_M2_result.png')
