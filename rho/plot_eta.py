#!/usr/bin/env python3
"""
plot_eta.py  --  Compute and plot RC factor eta vs Q^2

Usage:
    python plot_eta.py [-born born_events.lund] [-rc rc_events.lund]
"""

import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from lund_reader import read_lund, resolve_lund_dir

parser = argparse.ArgumentParser(description='Compute RC factor eta vs Q^2')
parser.add_argument('-version', default='', help='Model version (e.g. akushevich, harut) — sets lund_VERSION as default lund dir')
parser.add_argument('-born', default=None, help='Born LUND file')
parser.add_argument('-rc',   default=None, help='RC LUND file')
parser.add_argument('-nev',  default=0, type=int, help='Max events to read per file (default: 0 = all)')
args = parser.parse_args()
if args.born and args.rc:
    born_lund = args.born
    rc_lund   = args.rc
    _lund_dir = os.path.dirname(os.path.abspath(born_lund))
else:
    _lund_dir = resolve_lund_dir(f'lund_{args.version}') if args.version else 'lund'
    if not os.path.isdir(_lund_dir):
        import sys
        avail = [d for d in os.listdir('.') if d.startswith('lund_') and not d.endswith('.py') and (os.path.isdir(d) or os.path.isfile(d))]
        sys.exit(f"Error: lund directory '{_lund_dir}' not found.\n"
                 f"  Use -version to specify one of: {avail}\n"
                 f"  Example: python plot_eta.py -version akushevich")
    born_lund = args.born or f'{_lund_dir}/born_events.lund'
    rc_lund   = args.rc   or f'{_lund_dir}/rc_events.lund'
max_ev    = args.nev

def _stat_name(lund):
    return lund[:-5]+'_stat.dat' if lund.endswith('.lund') else lund+'_stat.dat'
born_stat = _stat_name(born_lund)
rc_stat   = _stat_name(rc_lund)

Mp = 0.938272

def read_sigma(fname):
    """Read total cross section (sigma_nb) from a gen_stat.dat file."""
    with open(fname) as f:
        for line in f:
            if 'sigma_nb' in line:
                return float(line.split('=')[1].split()[0])
    raise RuntimeError(f"sigma_nb not found in {fname}")

# Q2 bins covering the full generator range [1.0, 9.0] GeV²
Q2_bins = np.array([1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5, 9.0])
Q2_centers = 0.5*(Q2_bins[:-1] + Q2_bins[1:])

def read_Q2(fname, max_ev=0):
    """Read Q2 values from any LUND file using universal reader."""
    Q2s, _, _, _, _ = read_lund(fname, max_ev=max_ev)
    return Q2s

# Read both files
print(f"Input files:")
print(f"  LUND (Born) : {born_lund}")
print(f"  LUND (RC)   : {rc_lund}")
print(f"  Stat (Born) : {born_stat}")
print(f"  Stat (RC)   : {rc_stat}")
print()
print(f"Reading {born_lund} ...")
Q2_born = read_Q2(born_lund, max_ev)
print(f"  {len(Q2_born)} events, Q2: {Q2_born.min():.2f}-{Q2_born.max():.2f}")

print(f"Reading {rc_lund} ...")
Q2_rc = read_Q2(rc_lund, max_ev)
print(f"  {len(Q2_rc)} events, Q2: {Q2_rc.min():.2f}-{Q2_rc.max():.2f}")

# Read total cross sections from stat files for correct normalization.
# Each accepted event carries weight = sigma_total / N_ev, so:
#   eta(bin) = [N_rc(bin) * sigma_rc / N_ev_rc] / [N_born(bin) * sigma_born / N_ev_born]
sigma_born = read_sigma(born_stat)
sigma_rc   = read_sigma(rc_stat)
print(f"\n  sigma_born = {sigma_born:.4e} nb")
print(f"  sigma_rc   = {sigma_rc:.4e} nb")
print(f"  sigma_rc/sigma_born = {sigma_rc/sigma_born:.4f}  (expected eta)")

# Histogram both in same Q2 bins
n_born, _ = np.histogram(Q2_born, bins=Q2_bins)
n_rc,   _ = np.histogram(Q2_rc,   bins=Q2_bins)

norm_born = len(Q2_born)
norm_rc   = len(Q2_rc)

eta     = np.zeros(len(Q2_centers))
eta_err = np.zeros(len(Q2_centers))

for i in range(len(Q2_centers)):
    nb = n_born[i]
    nr = n_rc[i]
    if nb > 0 and nr > 0:
        # eta = (nr * sigma_rc / N_ev_rc) / (nb * sigma_born / N_ev_born)
        eta[i] = (nr * norm_born * sigma_rc) / (nb * norm_rc * sigma_born)
        # Poisson errors
        rel_err_born = 1./np.sqrt(nb)
        rel_err_rc   = 1./np.sqrt(nr)
        eta_err[i]   = eta[i] * np.sqrt(rel_err_born**2 + rel_err_rc**2)
    else:
        eta[i]     = np.nan
        eta_err[i] = np.nan

# Print table
print(f"\n{'Q2':>6} {'N_born':>8} {'N_rc':>8} {'eta':>8} {'err':>8}")
for i in range(len(Q2_centers)):
    print(f"{Q2_centers[i]:6.2f} {n_born[i]:8d} {n_rc[i]:8d} "
          f"{eta[i]:8.4f} {eta_err[i]:8.4f}")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))

# Our MC result
mask = np.isfinite(eta)
ax.errorbar(Q2_centers[mask], eta[mask], yerr=eta_err[mask],
            fmt='o', color='black', markerfacecolor='black',
            markersize=7, capsize=3, linewidth=1.5,
            label=r'MC generator ($v_{\rm max}=1.2$ GeV$^2$)')

# Load idiffrad result if available
if os.path.exists('result_vmax12.txt'):
    data = []
    with open('result_vmax12.txt') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 8:
                try:
                    data.append([float(parts[2]), float(parts[5])])
                except: pass
    if data:
        data = np.array(data)
        ax.plot(data[:,0], data[:,1], 's--', color='red',
                markersize=6, markerfacecolor='none',
                label=r'idiffrad exact ($v_{\rm max}=1.2$ GeV$^2$)')

ax.axhline(1.0, color='gray', linestyle=':', linewidth=1)
ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=13)
ax.set_ylabel(r'$\eta = \sigma_{\rm obs}/\sigma_{\rm Born}$', fontsize=13)
ax.set_xlim(0.8, 9.5)
ax.set_xscale('log')
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_ylim(0.75, 1.15)
ax.set_title(
    r'RC factor $\rho^0$ electroproduction, CLAS12 kinematics' + '\n' +
    r'$E_{\rm beam}=10.6$ GeV, $\sqrt{S}=4.46$ GeV, $v_{\rm max}=1.2$ GeV$^2$',
    fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, linestyle='--', alpha=0.4)
ax.tick_params(direction='in', which='both')

fig.tight_layout()
fig.savefig(f'{_lund_dir}/eta_mc_vs_exact.png', dpi=150)
print(f"\nSaved: {_lund_dir}/eta_mc_vs_exact.png")
os.system(f'open {_lund_dir}/eta_mc_vs_exact.png')
