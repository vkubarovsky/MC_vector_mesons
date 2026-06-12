#!/usr/bin/env python3
"""
compare_lund.py  --  Overlay kinematic distributions from two LUND files.
Works with any LUND format (DIFFRAD, gagrho, etc.).

Usage:
    python compare_lund.py -f1 file1.lund -f2 file2.lund [-l1 label1] [-l2 label2] [-o out.png]
"""
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from lund_reader import read_lund

parser = argparse.ArgumentParser()
parser.add_argument('-f1',  default='born_events.lund',  help='First LUND file')
parser.add_argument('-f2',  default='born_events.lund',  help='Second LUND file')
parser.add_argument('-l1',  default=None, help='Label for first file (default: filename)')
parser.add_argument('-l2',  default=None, help='Label for second file (default: filename)')
parser.add_argument('-o',   default=None, help='Output PNG (default: compare_lund.png)')
parser.add_argument('-nev', default=0, type=int, help='Max events per file (0=all)')
args = parser.parse_args()

label1 = args.l1 or os.path.basename(args.f1)
label2 = args.l2 or os.path.basename(args.f2)
outfile = args.o or 'compare_lund.png'

print(f"Reading {args.f1}  [{label1}]")
Q2a, xBa, ta, Wa, Epa = read_lund(args.f1, max_ev=args.nev)
print(f"  N={len(Q2a)}, Q2={Q2a.mean():.3f}+/-{Q2a.std():.3f}, xB={xBa.mean():.3f}, <t>={ta.mean():.3f}, W={Wa.mean():.3f}, E'={Epa.mean():.3f}")

print(f"Reading {args.f2}  [{label2}]")
Q2b, xBb, tb, Wb, Epb = read_lund(args.f2, max_ev=args.nev)
print(f"  N={len(Q2b)}, Q2={Q2b.mean():.3f}+/-{Q2b.std():.3f}, xB={xBb.mean():.3f}, <t>={tb.mean():.3f}, W={Wb.mean():.3f}, E'={Epb.mean():.3f}")

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(f'{label1}  vs  {label2}', fontsize=12, fontweight='bold')

def plot_comp(ax, xa, xb, bins, xrange, xlabel, log=True):
    kw = dict(bins=bins, range=xrange, density=True)
    ax.hist(xb, color='orange',    alpha=0.5, label=f'{label2}  mean={np.mean(xb):.3f}', **kw)
    ax.hist(xa, color='steelblue', alpha=0.5, label=f'{label1}  mean={np.mean(xa):.3f}', **kw)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel('Density (norm)', fontsize=10)
    if log: ax.set_yscale('log')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(direction='in', which='both')

plot_comp(axes[0,0], Q2a, Q2b, bins=40, xrange=(1, 6),
          xlabel=r'$Q^2$ (GeV$^2$)', log=False)
axes[0,0].set_title(r'$Q^2$ (linear scale)')

plot_comp(axes[0,1], Q2a, Q2b, bins=40, xrange=(1, 6),
          xlabel=r'$Q^2$ (GeV$^2$)', log=True)
axes[0,1].set_title(r'$Q^2$ (log scale)')

plot_comp(axes[0,2], xBa, xBb, bins=40, xrange=(0.05, 0.70),
          xlabel=r'$x_B$', log=True)
axes[0,2].set_title(r'$x_B$')

plot_comp(axes[1,0], ta,  tb,  bins=40, xrange=(-3.0, 0.0),
          xlabel=r'$t$ (GeV$^2$)', log=True)
axes[1,0].set_title(r'$t$ distribution')

plot_comp(axes[1,1], Wa,  Wb,  bins=40, xrange=(2.0, 5.5),
          xlabel=r'$W$ (GeV)', log=False)
axes[1,1].set_title(r'$W$ distribution')

plot_comp(axes[1,2], Epa, Epb, bins=40, xrange=(1.0, 10.5),
          xlabel=r"$E'$ (GeV)", log=False)
axes[1,2].set_title(r"Scattered lepton $E'$")

fig.tight_layout()
fig.savefig(outfile, dpi=150)
print(f"Saved: {outfile}")
os.system(f'open {outfile}')
