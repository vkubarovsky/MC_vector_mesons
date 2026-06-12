#!/usr/bin/env python3
"""
fit_params.py — Fit DIFFRAD J/psi parameters to match Mariana distributions.

Re-weights existing DIFFRAD Born events with trial σ_T parameters and
minimizes chi² against Mariana histograms in Q², W, -t, y, xB.

Usage:
    python fit_params.py -diffrad born_events.lund -mariana Mariana.lund
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import argparse, os
from lund_reader import read_lund_events

parser = argparse.ArgumentParser()
parser.add_argument('-diffrad', required=True, help='DIFFRAD Born LUND file')
parser.add_argument('-mariana', required=True, help='Mariana LUND file')
parser.add_argument('-nev', type=int, default=0, help='Max events per file (0=all)')
args = parser.parse_args()

Mp   = 0.938272
mJ   = 3.0969
mJsq = mJ**2
Wth  = Mp + mJ
Wth2 = Wth**2

def load_kin(fname, nev):
    print(f"Reading {fname} ...")
    evts = read_lund_events(fname, max_ev=nev)
    Q2 = np.array([e['Q2'] for e in evts])
    W  = np.array([e['W']  for e in evts])
    t  = np.array([e['t']  for e in evts])
    y  = np.array([e['y']  for e in evts])
    xB = np.array([e['xB'] for e in evts])
    print(f"  {len(evts)} events")
    return Q2, W, t, y, xB

dQ2, dW, dt, dy, dxB = load_kin(args.diffrad, args.nev)
mQ2, mW, mt, my, mxB = load_kin(args.mariana, args.nev)

# ── Current DIFFRAD parameters ───────────────────────────────────────────────
P0 = dict(alf1=400., alf2=1.0, alf3=0.32, nuT=3.0, mg2=1.0, cR=0.4)

def full_weight(Q2, W, t, y, p):
    W2 = W**2
    ok = W2 > Wth2
    sigT = np.zeros_like(Q2)
    sigT[ok] = (p['alf1']
                * (1. - Wth2/W2[ok])**p['alf2']
                * np.sqrt(W2[ok])**p['alf3']
                / (1. + Q2[ok]/mJsq)**p['nuT'])
    dsdt_T = sigT * 3. * p['mg2']**3 / (p['mg2'] - t)**4
    gamma2 = np.where(Q2 > 0, 4.*Mp**2 * (Q2/(W2 - Mp**2 + Q2))**2 / Q2, 0.)
    eps_n = 1. - y - 0.25*y**2*gamma2
    eps_d = 1. - y + 0.5*y**2 + 0.25*y**2*gamma2
    eps = np.where(eps_d > 0, np.clip(eps_n/eps_d, 0., 1.), 0.)
    return dsdt_T * (1. + eps * p['cR'] * Q2/mJsq)

sig0 = full_weight(dQ2, dW, dt, dy, P0)
sig0 = np.maximum(sig0, 1e-30)

# ── Histogram setup ──────────────────────────────────────────────────────────
bins_Q2 = np.linspace(0.001, 0.25, 25)
bins_W  = np.linspace(4.05, 4.56, 25)
bins_mt = np.linspace(0.0, 5.0, 25)
bins_y  = np.linspace(0.70, 1.0, 25)
bins_xB = np.linspace(0.0, 0.016, 25)

all_bins = [bins_Q2, bins_W, bins_mt, bins_y, bins_xB]
d_vars   = [dQ2, dW, -dt, dy, dxB]
m_vars   = [mQ2, mW, -mt, my, mxB]
labels   = [r'$Q^2$ (GeV$^2$)', r'$W$ (GeV)', r'$-t$ (GeV$^2$)', r'$y$', r'$x_B$']

mar_hists = []
for v, b in zip(m_vars, all_bins):
    h, _ = np.histogram(v, bins=b)
    h = h.astype(float)
    s = h.sum()
    if s > 0: h /= s
    mar_hists.append(h)

def chi2(params):
    alf1, alf2, alf3, nuT, mg2, cR = params
    p = dict(alf1=alf1, alf2=alf2, alf3=alf3, nuT=nuT, mg2=mg2, cR=cR)
    sig_new = full_weight(dQ2, dW, dt, dy, p)
    w = sig_new / sig0

    total = 0.
    for var, bins, h_mar in zip(d_vars, all_bins, mar_hists):
        h_dif, _ = np.histogram(var, bins=bins, weights=w)
        s = h_dif.sum()
        if s > 0: h_dif /= s
        ok = h_mar > 0
        if ok.sum() == 0: continue
        diff = h_dif[ok] - h_mar[ok]
        err2 = h_mar[ok] / len(mQ2) + 1e-6
        total += np.sum(diff**2 / err2)
    return total

# ── Fit using differential evolution (global, respects bounds) ───────────────
bounds = [(1., 1e7),     # alf1
          (0.1, 4.0),    # alf2
          (-2.0, 8.0),   # alf3
          (1.0, 5.0),    # nuT
          (0.3, 8.0),    # mg2
          (0.0, 1.5)]    # cR

print(f"\nInitial chi2 = {chi2([P0['alf1'], P0['alf2'], P0['alf3'], P0['nuT'], P0['mg2'], P0['cR']]):.1f}")
print("Running differential evolution (may take ~1 min) ...")

result = differential_evolution(chi2, bounds, seed=42, maxiter=500,
                                tol=1e-6, polish=True, disp=True)

pf = result.x
print(f"\n{'='*50}")
print(f"Best fit (chi2 = {result.fun:.1f}):")
print(f"  alf1 = {pf[0]:.1f}")
print(f"  alf2 = {pf[1]:.4f}")
print(f"  alf3 = {pf[2]:.4f}")
print(f"  nuT  = {pf[3]:.4f}")
print(f"  mg2  = {pf[4]:.4f}")
print(f"  cR   = {pf[5]:.4f}")
print(f"{'='*50}")

# ── Fortran snippet ──────────────────────────────────────────────────────────
print(f"\nFortran parameters for sigma_T_jpsi / sigma_L_jpsi:")
print(f"      parameter( alf1 = {pf[0]:.1f}d0  )")
print(f"      parameter( alf2 =   {pf[1]:.4f}d0 )")
print(f"      parameter( alf3 =   {pf[2]:.4f}d0)")
print(f"      parameter( nuT  =   {pf[3]:.4f}d0 )")
print(f"      parameter( mg2  =   {pf[4]:.4f}d0 )")
print(f"      parameter( cR   =   {pf[5]:.4f}d0 )")

# ── Plot ─────────────────────────────────────────────────────────────────────
p_best = dict(alf1=pf[0], alf2=pf[1], alf3=pf[2],
              nuT=pf[3], mg2=pf[4], cR=pf[5])
sig_best = full_weight(dQ2, dW, dt, dy, p_best)
w_best = sig_best / sig0

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle(r'DIFFRAD J/$\psi$ parameter fit to Mariana', fontsize=14, fontweight='bold')

for idx, (var_d, var_m, bins, lab) in enumerate(
        zip(d_vars, m_vars, all_bins, labels)):
    ax = axes.flat[idx]
    ax.hist(var_m, bins=bins, density=True,
            histtype='stepfilled', color='lightsalmon', edgecolor='darkred',
            alpha=0.7, label='Mariana')
    ax.hist(var_d, bins=bins, density=True,
            histtype='step', color='navy', lw=1.5, ls='--',
            label='DIFFRAD original')
    ax.hist(var_d, bins=bins, weights=w_best, density=True,
            histtype='step', color='blue', lw=2.5,
            label='DIFFRAD fitted')
    ax.set_xlabel(lab, fontsize=12)
    ax.set_ylabel('Normalized', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

ax = axes.flat[5]
ax.axis('off')
lines = [f"{'Parameter':<8} {'Original':>10} {'Fitted':>10}",
         f"{'─'*32}",
         f"{'α₁':<8} {P0['alf1']:>10.1f} {pf[0]:>10.1f}",
         f"{'α₂':<8} {P0['alf2']:>10.3f} {pf[1]:>10.4f}",
         f"{'α₃':<8} {P0['alf3']:>10.3f} {pf[2]:>10.4f}",
         f"{'νT':<8} {P0['nuT']:>10.3f} {pf[3]:>10.4f}",
         f"{'mg²':<8} {P0['mg2']:>10.3f} {pf[4]:>10.4f}",
         f"{'cR':<8} {P0['cR']:>10.3f} {pf[5]:>10.4f}",
         f"",
         f"χ² = {result.fun:.1f}"]
ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes, fontsize=13,
        fontfamily='monospace', va='top',
        bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', ec='gray'))

fig.tight_layout()
outfile = os.path.join(os.path.dirname(os.path.abspath(args.diffrad)), 'fit_params_result.png')
fig.savefig(outfile, dpi=150)
print(f"\nPlot saved: {outfile}")
os.system(f'open {outfile}')
