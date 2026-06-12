#!/usr/bin/env python3
"""
plot_eta2.py  --  Compute and plot RC factor eta vs Q^2 and vs W.

Generalized version of plot_eta.py: binning is derived from the data
range instead of hardcoded FT-range bins, so it works for any channel
(J/psi FT: Q2 ~ 0.01-0.12; phi CLAS12: Q2 ~ 0.4-10; etc.).
The original plot_eta.py is unchanged.

Usage:
    python3 plot_eta2.py -born born_events.lund -rc rc_events.lund [-nev N] [-o out.png]
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, argparse

parser = argparse.ArgumentParser(description='RC factor eta vs Q2 and W (auto-binned)')
parser.add_argument('-born', required=True, help='Born LUND file')
parser.add_argument('-rc',   required=True, help='RC LUND file')
parser.add_argument('-nev',  default=0, type=int, help='Max events per file (0 = all)')
parser.add_argument('-nbins', default=12, type=int, help='Number of bins per axis')
parser.add_argument('-o', default='eta_q2_w.png', help='Output PNG')
args = parser.parse_args()


def stat_name(lund):
    return lund[:-5] + '_stat.dat' if lund.endswith('.lund') else lund + '_stat.dat'


def read_sigma(fname):
    sig = err = None
    with open(fname) as f:
        for line in f:
            if 'sigma_nb' in line:
                sig = float(line.split('=')[1].split()[0])
            if 'sigma_err' in line:
                err = float(line.split('=')[1].split()[0])
    if sig is None:
        raise RuntimeError(f'sigma_nb not found in {fname}')
    return sig, (err or 0.0)


def read_headers(fname, max_ev=0):
    """W and Q2 from LUND headers (cols 13, 14; beam pid col 6 == 11)."""
    W2s, Q2s = [], []
    with open(fname) as f:
        for line in f:
            c = line.split()
            if len(c) == 15 and c[5] == '11':
                W2s.append(float(c[12]))
                Q2s.append(float(c[13]))
                if max_ev and len(W2s) >= max_ev:
                    break
    return np.sqrt(np.array(W2s)), np.array(Q2s)


def auto_bins(vals_b, vals_r, n):
    """Bin edges over the joint 0.5-99.5 percentile range; log spacing
    when the dynamic range exceeds one decade."""
    lo = min(np.percentile(vals_b, 0.5), np.percentile(vals_r, 0.5))
    hi = max(np.percentile(vals_b, 99.5), np.percentile(vals_r, 99.5))
    if lo > 0 and hi / lo > 10:
        return np.geomspace(lo, hi, n + 1), 'log'
    return np.linspace(lo, hi, n + 1), 'lin'


Wb, Q2b = read_headers(args.born, args.nev)
Wr, Q2r = read_headers(args.rc, args.nev)
sigB, sigB_e = read_sigma(stat_name(args.born))
sigR, sigR_e = read_sigma(stat_name(args.rc))
eta_tot = sigR / sigB
print(f'N_born={len(Wb)}  N_rc={len(Wr)}')
print(f'sigma_Born={sigB:.4e}  sigma_RC={sigR:.4e}  eta_tot={eta_tot:.4f}')

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, (vb, vr, label) in zip(axes, [(Wb, Wr, 'W [GeV]'),
                                      (Q2b, Q2r, r'$Q^2$ [GeV$^2$]')]):
    bins, scale = auto_bins(vb, vr, args.nbins)
    hB, _ = np.histogram(vb, bins=bins)
    hR, _ = np.histogram(vr, bins=bins)
    ctr = 0.5 * (bins[1:] + bins[:-1])
    dB = sigB * hB / len(vb)
    dR = sigR * hR / len(vr)
    m = hB > 20
    eta = np.where(m, dR / np.where(dB > 0, dB, 1), np.nan)
    err = eta * np.sqrt(1 / np.maximum(hB, 1) + 1 / np.maximum(hR, 1))
    ax.errorbar(ctr[m], eta[m], yerr=err[m], fmt='ko', ms=4)
    ax.axhline(eta_tot, color='r', ls='--',
               label=f'$\\eta_{{tot}}$={eta_tot:.3f}')
    ax.set_xlabel(label)
    ax.set_ylabel(r'$\eta$')
    ax.legend()
    if scale == 'log':
        ax.set_xscale('log')
    print(f'\n  {label} bins ({scale}):')
    for i in range(len(ctr)):
        if m[i]:
            print(f'   {bins[i]:9.4f}-{bins[i+1]:9.4f}  NB={hB[i]:7d}  '
                  f'NR={hR[i]:7d}  eta={eta[i]:.4f} +- {err[i]:.4f}')

fig.suptitle(os.path.basename(os.path.dirname(os.path.abspath(args.born))) or 'eta')
fig.tight_layout()
fig.savefig(args.o, dpi=130)
print(f'\nSaved: {args.o}')
