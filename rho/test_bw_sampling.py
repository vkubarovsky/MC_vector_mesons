#!/usr/bin/env python3
"""
test_bw_sampling.py  --  Reproduce sample_bw in Python and compare with
                         the theoretical running-width Breit-Wigner curve.

This tests the BW sampling algorithm from diffrad, WITHOUT any cross-section
accept/reject.  If the histogram matches the green BW curve, sample_bw is
correct.  Any remaining discrepancy in the LUND histogram is then due to
the M-dependent cross section in difflt (amv**3 factor etc.).

Usage:  python test_bw_sampling.py [-nev N] [-ammax AMMAX]
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nev',   default=1000000, type=int,
                    help='Number of samples (default: 1 000 000)')
parser.add_argument('-ammax', default=None, type=float,
                    help='Upper mass limit in GeV (default: M0+5*G0 = 1.526 GeV)')
args = parser.parse_args()

# ── BW parameters (from setcon, ivec=1) ──────────────────────────────────────
M0   = 0.77526   # rho nominal mass  (GeV)
G0   = 0.1502    # rho width         (GeV)
mpi  = 0.13957   # pion mass         (GeV)
ammin = 2.0 * mpi                           # threshold
ammax = args.ammax or M0 + 5.0 * G0        # upper limit (1.5258 GeV)

print(f"sample_bw simulation")
print(f"  M0    = {M0:.5f} GeV,  G0 = {G0:.5f} GeV,  mpi = {mpi:.5f} GeV")
print(f"  ammin = {ammin:.4f} GeV,  ammax = {ammax:.4f} GeV")
print(f"  N_try = {args.nev:,}")

rng  = np.random.default_rng(42)

p0     = np.sqrt(max(0.0, M0**2 / 4.0 - mpi**2))
thmin  = np.arctan((ammin**2 - M0**2) / (M0 * G0))
thmax  = np.arctan((ammax**2 - M0**2) / (M0 * G0))

# ── Simulate sample_bw ────────────────────────────────────────────────────────
samples = []
n_try   = 0
n_acc   = 0

while n_acc < args.nev:
    n_try += 1
    theta = thmin + rng.random() * (thmax - thmin)
    amv2  = M0**2 + M0 * G0 * np.tan(theta)
    amv   = np.sqrt(max(ammin**2, amv2))

    # running-width correction (ivec == 1, rho only)
    pm     = np.sqrt(max(0.0, amv**2 / 4.0 - mpi**2))
    gamrun = G0 * (pm / p0)**3 * (M0 / amv)
    bw_run   = 1.0 / ((amv2 - M0**2)**2 + M0**2 * gamrun**2)
    bw_const = 1.0 / ((amv2 - M0**2)**2 + M0**2 * G0**2)
    ratio    = bw_run / bw_const * (gamrun / G0)

    if rng.random() <= ratio:
        samples.append(amv)
        n_acc += 1

samples  = np.array(samples)
eff      = n_acc / n_try
print(f"  Accepted {n_acc:,} / {n_try:,}  (efficiency = {100*eff:.1f}%)")
print(f"  M range in sample: {samples.min():.4f} – {samples.max():.4f} GeV")
print(f"  Mean = {samples.mean():.4f} GeV,  Sigma = {samples.std():.4f} GeV")

# ── Theoretical BW curve (same formula) ──────────────────────────────────────
def bw_rho_dNdM(M):
    M2  = M**2
    pm  = np.sqrt(np.maximum(0.0, M2 / 4.0 - mpi**2))
    Gr  = G0 * (pm / p0)**3 * (M0 / M)
    return 2.0 * M * Gr / ((M2 - M0**2)**2 + M0**2 * Gr**2)

# ── Plot in histogram range (0.60, 0.95) and full range ──────────────────────
for (plot_range, nbins, tag) in [
        ((0.60, 0.95), 60,  '_zoom'),
        ((ammin, min(ammax, 1.30)), 80, '_full')]:
    Mmin, Mmax = plot_range
    bw_width   = (Mmax - Mmin) / nbins

    counts, edges = np.histogram(samples, bins=nbins, range=(Mmin, Mmax))
    N_in = counts.sum()

    M_line  = np.linspace(Mmin, Mmax, 800)
    bw_vals = bw_rho_dNdM(M_line)
    integral = np.trapz(bw_vals, M_line)
    scale    = N_in * bw_width / integral

    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(centers, counts, width=bw_width, color='lightskyblue',
           edgecolor='navy', alpha=0.8, label=f'sample_bw simulation ({n_acc:,} events)')
    ax.plot(M_line, scale * bw_vals, '-', color='forestgreen', lw=2.5,
            label=r'Theory: $2M\Gamma_{\rm run}/[(M^2-M_0^2)^2+M_0^2\Gamma_{\rm run}^2]$')
    ax.axvline(M0, color='red', lw=1.5, ls='--',
               label=f'$M_0$ = {M0:.5f} GeV')

    ax.set_xlabel(r'$M_{\pi^+\pi^-}$ (GeV)', fontsize=13)
    ax.set_ylabel('Events', fontsize=12)
    ax.set_title(
        r'BW sampling test: sample\_bw vs theory' + '\n' +
        rf'$M_0$={M0}, $\Gamma_0$={G0}, range [{Mmin:.3f}, {Mmax:.3f}] GeV',
        fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.tick_params(direction='in', which='both')

    ax.text(0.97, 0.95,
            f'N = {n_acc:,}\nMean = {samples[(samples>=Mmin)&(samples<=Mmax)].mean():.4f} GeV\n'
            f'Sigma = {samples[(samples>=Mmin)&(samples<=Mmax)].std():.4f} GeV',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow', ec='gray', alpha=0.9))

    fig.tight_layout()
    out = f'bw_sampling_test{tag}.png'
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")

import os
os.system('open bw_sampling_test_zoom.png bw_sampling_test_full.png')
