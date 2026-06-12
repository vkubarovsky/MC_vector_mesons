#!/usr/bin/env python3
"""
fit_Wdep.py — Fit α₂, α₃ in the W-dependence of σ_T to match Mariana.

DIFFRAD:  σ_T(W) ∝ (1 - W²_th/W²)^α₂ × (√W²)^α₃
Current:  α₂ = 1.0,  α₃ = 0.32

Re-weights DIFFRAD Born events and compares W histogram to Mariana.

Usage:
    python fit_Wdep.py -diffrad born_events.lund -mariana Mariana.lund
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brute
import argparse, os
from lund_reader import read_lund_events

parser = argparse.ArgumentParser()
parser.add_argument('-diffrad', required=True)
parser.add_argument('-mariana', required=True)
parser.add_argument('-nev', type=int, default=0)
args = parser.parse_args()

Mp   = 0.938272
mJ   = 3.0969
Wth  = Mp + mJ
Wth2 = Wth**2

def load_W(fname, nev):
    print(f"Reading {fname} ...")
    evts = read_lund_events(fname, max_ev=nev)
    W = np.array([e['W'] for e in evts])
    print(f"  {len(evts)} events,  W range: [{W.min():.3f}, {W.max():.3f}]")
    return W

dW = load_W(args.diffrad, args.nev)
mW = load_W(args.mariana, args.nev)

# ── Current parameters ───────────────────────────────────────────────────────
alf2_0, alf3_0 = 1.0, 0.32

def sigma_W(W, alf2, alf3):
    """W-dependent part of σ_T (without normalization)."""
    W2 = W**2
    thr = np.clip(1.0 - Wth2/W2, 0., None)
    return thr**alf2 * np.sqrt(W2)**alf3

# Reference weights for DIFFRAD events
sig0 = sigma_W(dW, alf2_0, alf3_0)
sig0 = np.maximum(sig0, 1e-30)

# ── Mariana W histogram ──────────────────────────────────────────────────────
bins_W = np.linspace(4.05, 4.56, 35)
h_mar, _ = np.histogram(mW, bins=bins_W)
h_mar = h_mar.astype(float)
norm_mar = h_mar.sum()
if norm_mar > 0:
    h_mar /= norm_mar

def chi2_W(params):
    alf2, alf3 = params
    if alf2 < 0.01:
        return 1e12
    sig_new = sigma_W(dW, alf2, alf3)
    w = sig_new / sig0
    h_dif, _ = np.histogram(dW, bins=bins_W, weights=w)
    s = h_dif.sum()
    if s > 0:
        h_dif /= s
    ok = h_mar > 0
    diff = h_dif[ok] - h_mar[ok]
    err2 = h_mar[ok] / norm_mar + 1e-6
    return np.sum(diff**2 / err2)

# ── Grid scan ────────────────────────────────────────────────────────────────
a2_grid = np.linspace(0.1, 4.0, 80)
a3_grid = np.linspace(-2.0, 4.0, 80)
chi2_map = np.zeros((len(a3_grid), len(a2_grid)))
for i, a3 in enumerate(a3_grid):
    for j, a2 in enumerate(a2_grid):
        chi2_map[i, j] = chi2_W([a2, a3])

# ── Optimize from several starting points ────────────────────────────────────
best = None
for a2_start in [0.3, 0.5, 1.0, 1.5, 2.0, 3.0]:
    for a3_start in [-1.0, 0.0, 0.32, 1.0, 2.0, 3.0]:
        res = minimize(chi2_W, [a2_start, a3_start], method='Nelder-Mead',
                       options=dict(xatol=1e-5, fatol=1e-8, maxiter=20000))
        if best is None or res.fun < best.fun:
            best = res

alf2_f, alf3_f = best.x
print(f"\nBest fit:  α₂ = {alf2_f:.4f},  α₃ = {alf3_f:.4f}")
print(f"  χ²(W) = {best.fun:.2f}")

# Also show some fixed-α₃ fits
print(f"\n{'α₂':>8s} {'α₃':>8s} {'χ²(W)':>10s}")
print("-" * 30)
configs = []
for a3_fix in [-1.0, -0.5, 0.0, 0.32, 0.5, 1.0, 2.0]:
    res_a2 = minimize(lambda x: chi2_W([x[0], a3_fix]), [1.0],
                      method='Nelder-Mead', options=dict(xatol=1e-5))
    c = chi2_W([res_a2.x[0], a3_fix])
    configs.append((res_a2.x[0], a3_fix, c))
    print(f"{res_a2.x[0]:8.4f} {a3_fix:8.2f} {c:10.2f}")
configs.append((alf2_f, alf3_f, best.fun))
print(f"{alf2_f:8.4f} {alf3_f:8.4f} {best.fun:10.2f}  ← free fit")
configs.append((alf2_0, alf3_0, chi2_W([alf2_0, alf3_0])))
print(f"{alf2_0:8.4f} {alf3_0:8.4f} {chi2_W([alf2_0, alf3_0]):10.2f}  ← current")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle(r'Fitting $W$-dependence: $(1-W_\mathrm{th}^2/W^2)^{\alpha_2}\,W^{\alpha_3}$',
             fontsize=14, fontweight='bold')

bc = 0.5*(bins_W[:-1] + bins_W[1:])

# (a) W distributions
ax = axes[0]
ax.hist(mW, bins=bins_W, density=True, histtype='stepfilled',
        color='lightsalmon', edgecolor='darkred', alpha=0.7, label='Mariana')
ax.hist(dW, bins=bins_W, density=True, histtype='step',
        color='navy', lw=1.5, ls='--', label=f'DIFFRAD orig (α₂={alf2_0}, α₃={alf3_0})')

# Re-weighted best fit
w_best = sigma_W(dW, alf2_f, alf3_f) / sig0
ax.hist(dW, bins=bins_W, weights=w_best, density=True, histtype='step',
        color='blue', lw=2.5, label=f'DIFFRAD fit (α₂={alf2_f:.2f}, α₃={alf3_f:.2f})')

# A few other interesting fits
nice_colors = ['green', 'orange', 'purple']
for idx, a3_show in enumerate([0.0, 1.0, 2.0]):
    # find best α₂ for this α₃
    for a2v, a3v, _ in configs:
        if abs(a3v - a3_show) < 0.01:
            w_c = sigma_W(dW, a2v, a3_show) / sig0
            ax.hist(dW, bins=bins_W, weights=w_c, density=True, histtype='step',
                    color=nice_colors[idx], lw=1.5,
                    label=f'α₂={a2v:.2f}, α₃={a3_show:.1f}')
            break

ax.set_xlabel(r'$W$ (GeV)', fontsize=12)
ax.set_ylabel('Normalized', fontsize=10)
ax.set_title('$W$ distribution comparison')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

# (b) σ_T(W) shape comparison
ax = axes[1]
Wplot = np.linspace(Wth + 0.001, 4.56, 200)
# Normalize all to max=1
for label, a2, a3, col, ls, lw in [
        (f'Current ({alf2_0}, {alf3_0})', alf2_0, alf3_0, 'navy', '--', 1.5),
        (f'Best fit ({alf2_f:.2f}, {alf3_f:.2f})', alf2_f, alf3_f, 'blue', '-', 2.5),
        (f'α₃=0', configs[2][0], 0.0, 'green', '-', 1.5),
        (f'α₃=1', configs[5][0], 1.0, 'orange', '-', 1.5),
        (f'α₃=2', configs[6][0], 2.0, 'purple', '-', 1.5)]:
    s = sigma_W(Wplot, a2, a3)
    s /= s.max()
    ax.plot(Wplot, s, color=col, ls=ls, lw=lw, label=label)
ax.set_xlabel(r'$W$ (GeV)', fontsize=12)
ax.set_ylabel(r'$\sigma_T(W)$ / max', fontsize=12)
ax.set_title(r'$\sigma_T(W)$ shape')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (c) χ² landscape
ax = axes[2]
vmin = np.log10(max(best.fun, 0.01))
vmax = np.log10(min(chi2_map.max(), 1e6))
cs = ax.contourf(a2_grid, a3_grid, np.log10(np.clip(chi2_map, 0.01, 1e6)),
                 levels=30, cmap='viridis', vmin=vmin, vmax=vmax)
plt.colorbar(cs, ax=ax, label=r'log$_{10}(\chi^2)$')
ax.plot(alf2_f, alf3_f, 'r*', ms=15, label=f'Best: ({alf2_f:.2f}, {alf3_f:.2f})')
ax.plot(alf2_0, alf3_0, 'wo', ms=10, mew=2, label=f'Current: ({alf2_0}, {alf3_0})')
ax.set_xlabel(r'$\alpha_2$', fontsize=12)
ax.set_ylabel(r'$\alpha_3$', fontsize=12)
ax.set_title(r'$\chi^2$ landscape ($W$ only)')
ax.legend(fontsize=9)

fig.tight_layout()
fig.savefig('fit_Wdep_result.png', dpi=150)
print(f"\nPlot saved: fit_Wdep_result.png")
os.system('open fit_Wdep_result.png')
