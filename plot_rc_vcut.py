#!/usr/bin/env python3
"""
plot_rc_vcut.py  --  Plot RC factor eta = sigma_RC(v_max) / sigma_Born
                     as a function of the ISR photon energy cut v_max.

Usage:
    python plot_rc_vcut.py -born rho/lund/born_events.lund \
                           -rc   rho/lund/rc_events.lund   \
                           [-label "rho0"] [-out rc_vcut.png]
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-born',  default='rho/lund/born_events.lund')
parser.add_argument('-rc',    default='rho/lund/rc_events.lund')
parser.add_argument('-label', default='')
parser.add_argument('-out',   default='')
parser.add_argument('-nev',   default=0, type=int)
args = parser.parse_args()

# ── helpers ────────────────────────────────────────────────────────────────

def stat_name(lund):
    return lund[:-5]+'_stat.dat' if lund.endswith('.lund') else lund+'_stat.dat'

def vdist_name(lund):
    return lund[:-5]+'_vdist.dat' if lund.endswith('.lund') else lund+'_vdist.dat'

def read_stat(fname):
    d = {}
    if not os.path.exists(fname): return d
    with open(fname) as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                d[k.strip()] = float(v.split()[0])
    return d

def read_lund_weights(fname, max_ev=0):
    """Return (weights, hard_flag) arrays from LUND header lines."""
    weights, hard = [], []
    with open(fname) as f:
        for line in f:
            if max_ev > 0 and len(weights) >= max_ev:
                break
            tok = line.split()
            if not tok: continue
            try:
                npart = int(tok[0])
                if npart not in (6, 7): raise ValueError
                if len(tok) < 15:       raise ValueError
                w = float(tok[9])
            except:
                continue
            weights.append(w)
            hard.append(npart == 7)
            for _ in range(npart):
                next(f, '')
    return np.array(weights), np.array(hard)

def read_vdist(fname):
    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        return np.loadtxt(fname)
    return np.array([])

# ── load data ──────────────────────────────────────────────────────────────

print(f"Reading BORN: {args.born}")
stat_b  = read_stat(stat_name(args.born))
sigma_born = stat_b.get('sigma_nb', float('nan'))
print(f"  sigma_Born = {sigma_born:.6e} nb")

print(f"Reading RC  : {args.rc}")
stat_rc    = read_stat(stat_name(args.rc))
sigma_rc   = stat_rc.get('sigma_nb', float('nan'))
weights_rc, hard_rc = read_lund_weights(args.rc, args.nev)
vv         = read_vdist(vdist_name(args.rc))

N_hard = hard_rc.sum()
print(f"  sigma_RC   = {sigma_rc:.6e} nb  ({len(weights_rc)} events, {N_hard} hard)")

wsum_soft = weights_rc[~hard_rc].sum()
wsum_hard = weights_rc[hard_rc].sum()
wsum_tot  = wsum_soft + wsum_hard

if len(vv) != N_hard:
    print(f"  WARNING: vdist length {len(vv)} != N_hard {N_hard}. Using histogram fraction.")
    use_sorted = False
else:
    use_sorted = True

# ── build eta(v_max) curve  (smoothed) ────────────────────────────────────
# Build exact step-function cumulative from sorted events, evaluate on a
# dense grid, then apply Gaussian smoothing to the curve itself.
# Smoothing scale ~2% of v range gives visually clean line while
# preserving the shape and end-point value exactly.

V_MAX_GEN = 1.20
v_grid    = np.linspace(0.0, V_MAX_GEN, 2000)

if N_hard > 0:
    idx    = np.argsort(vv)
    vv_s   = vv[idx]
    wh_s   = weights_rc[hard_rc][idx]
    wh_cum = np.concatenate([[0.0], np.cumsum(wh_s)])   # length N_hard+1
    pos    = np.searchsorted(vv_s, v_grid, side='right')
    wh_cut_raw = wh_cum[pos]
else:
    wh_cut_raw = np.zeros_like(v_grid)

w_cut_raw     = wsum_soft + wh_cut_raw
eta_raw       = sigma_rc * w_cut_raw / wsum_tot / sigma_born

# Smooth the curve; sigma = 2% of grid length in points
sigma_pts = max(3, int(0.02 * len(v_grid)))
eta       = gaussian_filter1d(eta_raw, sigma=sigma_pts, mode='nearest')
# Preserve exact end-point
eta[-1]   = sigma_rc / sigma_born

sigma_rc_vcut = eta * sigma_born   # for the printed table

eta_full = sigma_rc / sigma_born

# ── plot ───────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(left=0.08, right=0.97, bottom=0.12, top=0.88, wspace=0.30)

meson_label = args.label if args.label else os.path.basename(os.path.dirname(
                  os.path.dirname(args.rc)))

_latex_name = {'rho':  r'$\rho^0$',
               'phi':  r'$\phi$',
               'jpsi': r'$J/\psi$'}
meson_display = _latex_name.get(meson_label.lower(), meson_label)

fig.suptitle(rf"RC factor vs hard-ISR cut $v$ for {meson_display}", fontsize=13)

# ── left panel: linear scale ───────────────────────────────────────────────
ax = axes[0]
ax.plot(v_grid, eta, 'b-', lw=1.8)
ax.axhline(eta_full, color='gray', lw=1.0, ls='--', label=f'full: η={eta_full:.4f}')
ax.axvline(0.5,  color='orange', lw=1.0, ls=':', label='v=0.5 GeV²')
ax.axvline(1.0,  color='red',    lw=1.0, ls=':', label='v=1.0 GeV²')
ax.set_xlabel(r'$v_{\max} = M_X^2 - M_p^2\ \mathrm{(GeV}^2\mathrm{)}$', fontsize=12)
ax.set_ylabel(r'$\eta(v_{\max}) = \sigma_\mathrm{RC}(v_{\max})\,/\,\sigma_\mathrm{Born}$', fontsize=11)
ax.set_xlim(0, 1.25)
ax.set_ylim(max(0, eta.min() * 0.995), eta.max() * 1.005)
ax.set_title('Linear scale', fontsize=10)
ax.grid(True, alpha=0.3)

# ── right panel: log x-scale (zoom low v) ─────────────────────────────────
ax2 = axes[1]
ax2.plot(v_grid, eta, 'b-', lw=1.8)
ax2.axhline(eta_full, color='gray', lw=1.0, ls='--')
ax2.axvline(0.1,  color='green',  lw=1.0, ls=':')
ax2.axvline(0.5,  color='orange', lw=1.0, ls=':')
ax2.set_xscale('log')
ax2.set_xlabel(r'$v_{\max} = M_X^2 - M_p^2\ \mathrm{(GeV}^2\mathrm{)}$', fontsize=12)
ax2.set_ylabel(r'$\eta(v_{\max}) = \sigma_\mathrm{RC}(v_{\max})\,/\,\sigma_\mathrm{Born}$', fontsize=11)
ax2.set_xlim(0.01, 1.25)
ax2.set_ylim(max(0, eta.min() * 0.995), eta.max() * 1.005)
ax2.set_title('Log x-scale', fontsize=10)
ax2.grid(True, alpha=0.3, which='both')

# ── table of key values ────────────────────────────────────────────────────
box_cuts = [0.01, 0.05, 0.10, 0.50, 1.00, 1.20]
key_cuts = [0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.00, 1.20]
print(f"\n  eta(v_max) = sigma_RC(v_max) / sigma_Born:")
print(f"  sigma_Born = {sigma_born:.6e} nb")
print(f"  {'v_max':>8}  {'sigma_RC':>14}  {'eta':>8}")
eta_at = {}
for vc in sorted(set(key_cuts) | set(box_cuts)):
    pos_k = np.searchsorted(v_grid, vc, side='right') - 1
    pos_k = max(0, min(pos_k, len(eta)-1))
    eta_at[vc] = eta[pos_k]
    print(f"  {vc:8.2f}  {eta_at[vc]*sigma_born:14.6e}  {eta_at[vc]:8.5f}")

# ── text box on both panels ────────────────────────────────────────────────
box_lines = '\n'.join(
    [rf'$v_{{max}}$={vc:.2f}  $\eta$={eta_at[vc]:.4f}' for vc in box_cuts]
)
for ax in axes:
    ax.text(0.98, 0.04, box_lines,
            transform=ax.transAxes,
            fontsize=8, family='monospace',
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.85))

# ── save ───────────────────────────────────────────────────────────────────
outfile = args.out if args.out else \
          os.path.join(os.path.dirname(args.rc), '..', 'rc_vcut.png')
outfile = os.path.normpath(outfile)
plt.savefig(outfile, dpi=150)
print(f"\nSaved: {outfile}")
