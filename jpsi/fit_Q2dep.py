#!/usr/bin/env python3
"""
fit_Q2dep.py — Fit νT in the Q² dependence of σ_T to match Mariana.

DIFFRAD:  σ_T(Q²) ∝ 1/(1 + Q²/m²_J/ψ)^νT
Current:  νT = 3.0,  m_J/ψ = 3.0969 GeV

Re-weights DIFFRAD Born events and compares Q² histogram to Mariana.

Usage:
    python fit_Q2dep.py -diffrad born_events.lund -mariana Mariana.lund
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
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

nuT_0 = 3.0

def sigma_Q2(Q2, nuT):
    return 1.0 / (1.0 + Q2/mJsq)**nuT

sig0 = sigma_Q2(dQ2, nuT_0)

# ── Mariana Q² histogram ─────────────────────────────────────────────────────
bins_Q2 = np.linspace(0.001, 0.25, 35)
h_mar, _ = np.histogram(mQ2, bins=bins_Q2)
h_mar = h_mar.astype(float)
norm_mar = h_mar.sum()
if norm_mar > 0:
    h_mar /= norm_mar

def chi2_Q2(nuT):
    if nuT < 0.0:
        return 1e12
    sig_new = sigma_Q2(dQ2, nuT)
    w = sig_new / sig0
    h_dif, _ = np.histogram(dQ2, bins=bins_Q2, weights=w)
    s = h_dif.sum()
    if s > 0:
        h_dif /= s
    ok = h_mar > 0
    diff = h_dif[ok] - h_mar[ok]
    err2 = h_mar[ok] / norm_mar + 1e-6
    return np.sum(diff**2 / err2)

# ── Scan νT ──────────────────────────────────────────────────────────────────
nuT_scan = np.linspace(0.0, 20.0, 200)
chi2_scan = np.array([chi2_Q2(n) for n in nuT_scan])

# ── Minimize ─────────────────────────────────────────────────────────────────
res = minimize_scalar(chi2_Q2, bounds=(0.0, 50.0), method='bounded')
nuT_fit = res.x

print(f"\nBest fit:  νT = {nuT_fit:.2f}")
print(f"  χ²(Q²) = {res.fun:.2f}")
print(f"  Current νT=3:  χ²(Q²) = {chi2_Q2(3.0):.2f}")

# Check: how much does 1/(1+Q²/mJ²)^νT actually vary over 0.001-0.25?
print(f"\n  1/(1+Q²/mJ²)^νT at Q²=0.001:  νT=3 → {sigma_Q2(0.001, 3.0):.6f}"
      f"   νT={nuT_fit:.1f} → {sigma_Q2(0.001, nuT_fit):.6f}")
print(f"  1/(1+Q²/mJ²)^νT at Q²=0.25:   νT=3 → {sigma_Q2(0.25, 3.0):.6f}"
      f"   νT={nuT_fit:.1f} → {sigma_Q2(0.25, nuT_fit):.6f}")
print(f"  Ratio max/min:  νT=3 → {sigma_Q2(0.001, 3.0)/sigma_Q2(0.25, 3.0):.4f}"
      f"   νT={nuT_fit:.1f} → {sigma_Q2(0.001, nuT_fit)/sigma_Q2(0.25, nuT_fit):.4f}")

# Table of chi² vs νT
print(f"\n{'νT':>6s} {'χ²(Q²)':>10s} {'ratio(Q²min/Q²max)':>20s}")
print("-" * 40)
for nuT_try in [0, 1, 2, 3, 4, 5, 10, 15, 20, nuT_fit]:
    c = chi2_Q2(nuT_try)
    r = sigma_Q2(0.001, nuT_try) / sigma_Q2(0.25, nuT_try)
    tag = "  ← best" if abs(nuT_try - nuT_fit) < 0.01 else ""
    tag = "  ← current" if abs(nuT_try - 3.0) < 0.01 else tag
    print(f"{nuT_try:6.1f} {c:10.2f} {r:20.4f}{tag}")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle(r'Fitting $Q^2$ dependence: $1/(1+Q^2/m_{J/\psi}^2)^{\nu_T}$',
             fontsize=14, fontweight='bold')

# (a) Q² distributions
ax = axes[0]
ax.hist(mQ2, bins=bins_Q2, density=True, histtype='stepfilled',
        color='lightsalmon', edgecolor='darkred', alpha=0.7, label='Mariana')
ax.hist(dQ2, bins=bins_Q2, density=True, histtype='step',
        color='navy', lw=1.5, ls='--', label=f'DIFFRAD orig (νT={nuT_0})')

w_fit = sigma_Q2(dQ2, nuT_fit) / sig0
ax.hist(dQ2, bins=bins_Q2, weights=w_fit, density=True, histtype='step',
        color='blue', lw=2.5, label=f'DIFFRAD fit (νT={nuT_fit:.1f})')

for nuT_show, col in [(5, 'green'), (10, 'orange'), (20, 'purple')]:
    w_s = sigma_Q2(dQ2, nuT_show) / sig0
    ax.hist(dQ2, bins=bins_Q2, weights=w_s, density=True, histtype='step',
            color=col, lw=1.5, label=f'νT={nuT_show}')

ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel('Normalized', fontsize=10)
ax.set_title(r'$Q^2$ distribution')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (a') same in log
ax = axes[1]
ax.hist(mQ2, bins=bins_Q2, density=True, histtype='stepfilled',
        color='lightsalmon', edgecolor='darkred', alpha=0.7, label='Mariana')
ax.hist(dQ2, bins=bins_Q2, density=True, histtype='step',
        color='navy', lw=1.5, ls='--', label=f'DIFFRAD orig (νT={nuT_0})')
ax.hist(dQ2, bins=bins_Q2, weights=w_fit, density=True, histtype='step',
        color='blue', lw=2.5, label=f'DIFFRAD fit (νT={nuT_fit:.1f})')
for nuT_show, col in [(5, 'green'), (10, 'orange'), (20, 'purple')]:
    w_s = sigma_Q2(dQ2, nuT_show) / sig0
    ax.hist(dQ2, bins=bins_Q2, weights=w_s, density=True, histtype='step',
            color=col, lw=1.5, label=f'νT={nuT_show}')
ax.set_yscale('log')
ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=12)
ax.set_ylabel('Normalized', fontsize=10)
ax.set_title(r'$Q^2$ distribution (log)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (b) χ² vs νT
ax = axes[2]
ax.plot(nuT_scan, chi2_scan, 'b-', lw=2)
ax.axvline(nuT_fit, color='red', ls='--', lw=1.5, label=f'Best: νT={nuT_fit:.1f}')
ax.axvline(nuT_0, color='gray', ls=':', lw=1.5, label=f'Current: νT={nuT_0}')
ax.set_xlabel(r'$\nu_T$', fontsize=12)
ax.set_ylabel(r'$\chi^2(Q^2)$', fontsize=12)
ax.set_title(r'$\chi^2$ vs $\nu_T$')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('fit_Q2dep_result.png', dpi=150)
print(f"\nPlot saved: fit_Q2dep_result.png")
os.system('open fit_Q2dep_result.png')
