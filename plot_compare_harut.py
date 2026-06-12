#!/usr/bin/env python3
"""
plot_compare_harut.py  --  plot_compare-style plots for Harut's external LUND format.

External LUND particle order (7 particles per event):
  1 = beam electron
  2 = target proton
  3 = scattered electron  (k2)
  4 = recoil proton       (pp)
  5 = rho0 / vector meson (SKIP)
  6 = pi-                 (pim)
  7 = pi+                 (pip)

Header line: only tok[0]=npart is used; all kinematic fields ignored.
Particle line: idx charge type PDG parent daughter px py pz E mass xv yv zv

Usage:
    python plot_compare_harut.py -lund file.lund [-nev N] [-out prefix]
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-lund',   required=True,  help='Input LUND file')
parser.add_argument('-nev',    default=0, type=int)
parser.add_argument('-label',  default='')
parser.add_argument('-out',    default='',    help='Output file prefix')
parser.add_argument('-q2min',  default=0.0, type=float, help='Q2 minimum cut (GeV^2)')
parser.add_argument('-q2max',  default=0.0, type=float, help='Q2 maximum cut (GeV^2)')
parser.add_argument('-wmin',   default=0.0, type=float, help='W minimum cut (GeV)')
args = parser.parse_args()

Mp   = 0.938272
Mpi  = 0.139570

# ── LUND reader ────────────────────────────────────────────────────────────

def read_lund(fname, max_ev=0):
    """Read external LUND file. Returns list of dicts with 4-vectors."""
    events = []
    with open(fname) as f:
        while True:
            if max_ev > 0 and len(events) >= max_ev:
                break
            hdr = f.readline()
            if not hdr:
                break
            tok = hdr.split()
            if not tok:
                continue
            try:
                npart = int(tok[0])
                if npart < 7:
                    raise ValueError
            except:
                continue
            parts = []
            for _ in range(npart):
                line = f.readline()
                if not line:
                    break
                t = line.split()
                if len(t) < 11:
                    continue
                # (E, px, py, pz)
                parts.append(np.array([float(t[9]), float(t[6]),
                                        float(t[7]), float(t[8])]))
            if len(parts) < 7:
                continue
            events.append(parts)
    return events

def m2v(p4):
    return p4[:,0]**2 - p4[:,1]**2 - p4[:,2]**2 - p4[:,3]**2

def pkin(p4):
    pmag  = np.sqrt(p4[:,1]**2 + p4[:,2]**2 + p4[:,3]**2)
    theta = np.degrees(np.arctan2(np.sqrt(p4[:,1]**2 + p4[:,2]**2), p4[:,3]))
    phi   = np.degrees(np.arctan2(p4[:,2], p4[:,1]))
    return pmag, theta, phi

# ── extract kinematics ─────────────────────────────────────────────────────

def extract(events):
    k1  = np.array([ev[0] for ev in events])   # beam e-
    # ev[1] = target proton (not needed; use nominal)
    k2  = np.array([ev[2] for ev in events])   # scattered e-
    pp  = np.array([ev[3] for ev in events])   # recoil proton
    # ev[4] = rho  (skip)
    pim = np.array([ev[5] for ev in events])   # pi-
    pip = np.array([ev[6] for ev in events])   # pi+

    ptar = np.zeros_like(k1)
    ptar[:,0] = Mp                              # target at rest

    q4   = k1 - k2                             # virtual photon
    Q2   = -(q4[:,0]**2 - q4[:,1]**2 - q4[:,2]**2 - q4[:,3]**2)
    nu   = k1[:,0] - k2[:,0]
    y    = nu / k1[:,0]

    W2   = Mp**2 + 2*Mp*nu - Q2
    W    = np.sqrt(np.maximum(0, W2))

    ph   = pip + pim
    Mrho = np.sqrt(np.maximum(0, m2v(ph)))

    # t = (pp - ptar)^2
    dp   = pp - ptar
    t    = dp[:,0]**2 - dp[:,1]**2 - dp[:,2]**2 - dp[:,3]**2

    pT2  = ph[:,1]**2 + ph[:,2]**2
    xB   = Q2 / (2*Mp*nu)

    # missing mass: X = k1 + ptar - k2 - pip - pim - pp
    miss = k1 + ptar - k2 - pip - pim - pp
    Mmiss2 = m2v(miss)
    Mmiss  = np.sign(Mmiss2)*np.sqrt(np.abs(Mmiss2))
    E_miss = miss[:,0]

    return dict(Q2=Q2, y=y, t=t, W=W, Mrho=Mrho, Mmiss=Mmiss,
                E_miss=E_miss, nu=nu, xB=xB, pT2=pT2,
                k2_kin=pkin(k2), pp_kin=pkin(pp),
                pip_kin=pkin(pip), pim_kin=pkin(pim))

# ── annotation helper ──────────────────────────────────────────────────────

def annotate(ax, arr, unit='', loc='right'):
    x = 0.97 if loc == 'right' else 0.03
    ha = 'right' if loc == 'right' else 'left'
    txt = (f"N     = {len(arr)}\n"
           f"Mean = {np.mean(arr):.4f}{unit}\n"
           f"Sigma = {np.std(arr):.4f}{unit}")
    ax.text(x, 0.97, txt, transform=ax.transAxes,
            fontsize=7, va='top', ha=ha,
            bbox=dict(boxstyle='round', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.8))

# ── validation plot (3×4 grid) ─────────────────────────────────────────────

def make_plot(d, title, outfile):
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.07,
                        top=0.93, hspace=0.45, wspace=0.35)
    fig.suptitle(title, fontsize=13)

    kw = dict(bins=60, histtype='stepfilled', color='lightskyblue', edgecolor='navy')

    # row 0
    ax = axes[0,0]; ax.hist(d['Q2'], **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$Q^2$ (GeV$^2$)'); ax.set_ylabel('Events')
    ax.set_title(r'$Q^2$'); annotate(ax, d['Q2'])

    ax = axes[0,1]; ax.hist(d['xB'], **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$x_B$'); ax.set_ylabel('Events')
    ax.set_title(r'$x_B$'); annotate(ax, d['xB'])

    ax = axes[0,2]; ax.hist(-d['t'], **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$-t$ (GeV$^2$)'); ax.set_ylabel('Events')
    ax.set_title(r'$-t$'); annotate(ax, -d['t'])

    ax = axes[0,3]; ax.hist(d['y'], **kw)
    ax.set_xlabel(r'$y$'); ax.set_ylabel('Events')
    ax.set_title(r'$y$'); annotate(ax, d['y'], loc='left')

    # row 1
    ax = axes[1,0]; ax.hist(d['W'], **kw)
    ax.set_xlabel(r'$W$ (GeV)'); ax.set_ylabel('Events')
    ax.set_title(r'$W$'); annotate(ax, d['W'])

    ax = axes[1,1]
    mm = d['Mmiss']
    ax.hist(np.abs(mm), **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$|M_\mathrm{miss}|$ (GeV)'); ax.set_ylabel('Events')
    ax.set_title('Missing mass'); annotate(ax, mm)

    ax = axes[1,2]; ax.hist(d['pT2'], **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$p_T^2$ (GeV$^2$)'); ax.set_ylabel('Events')
    ax.set_title(r'Meson $p_T^2$'); annotate(ax, d['pT2'])

    ax = axes[1,3]; ax.hist(d['Mrho'], bins=80, range=(0.60, 0.95),
                            histtype='stepfilled', color='lightskyblue', edgecolor='navy')
    ax.set_xlabel(r'$M_{\pi^+\pi^-}$ (GeV)'); ax.set_ylabel('Events')
    ax.set_title(r'$\rho^0$ mass'); annotate(ax, d['Mrho'])

    # row 2 — placeholders for v/Eγ (not available), use xB/nu/E_miss/nu
    ax = axes[2,0]; ax.axis('off')   # no v dist in external format

    ax = axes[2,1]; ax.axis('off')   # no ISR photon

    ax = axes[2,2]; ax.hist(np.abs(d['E_miss']+1e-9), **kw)
    ax.set_yscale('log')
    ax.set_xlabel(r'$|E_\mathrm{miss}|$ (GeV)'); ax.set_ylabel('Events')
    ax.set_title(r'$E_\mathrm{miss}$'); annotate(ax, d['E_miss'])

    ax = axes[2,3]; ax.hist(d['nu'], **kw)
    ax.set_xlabel(r'$\nu$ (GeV)'); ax.set_ylabel('Events')
    ax.set_title(r'$\nu$'); annotate(ax, d['nu'], loc='left')

    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f'Saved: {outfile}')

# ── kinematics plot (4×3 grid, Born only) ─────────────────────────────────

def make_kinematics_plot(d, title, outfile,
                         plabels=(r"$e'$", r'$p$', r'$\pi^+$', r'$\pi^-$'),
                         theta_maxs=(20, 90, 40, 40)):
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.05,
                        top=0.93, hspace=0.45, wspace=0.35)
    fig.suptitle(title, fontsize=13)

    kin_keys = ['k2_kin', 'pp_kin', 'pip_kin', 'pim_kin']
    pmaxs    = [9, 1.5, 6, 6]

    for row, (key, pname, thmax, pmax) in enumerate(
            zip(kin_keys, plabels, theta_maxs, pmaxs)):
        pmag, theta, phi = d[key]

        # col 0: |p| histogram
        ax = axes[row, 0]
        ax.hist(pmag, bins=60, range=(0, pmax),
                histtype='stepfilled', color='lightskyblue', edgecolor='navy')
        ax.set_xlabel(r'$|p|$ (GeV/c)', fontsize=10)
        ax.set_ylabel('Events', fontsize=9)
        ax.set_title(f'{pname} momentum', fontsize=10)
        annotate(ax, pmag)

        # col 1: 2D theta vs |p|
        ax = axes[row, 1]
        h = ax.hist2d(pmag, theta, bins=60, cmap='plasma',
                      range=[[0, pmax], [0, thmax]])
        fig.colorbar(h[3], ax=ax)
        ax.set_xlabel(r'$|p|$ (GeV/c)', fontsize=10)
        ax.set_ylabel(r'$\theta$ (deg)', fontsize=10)
        ax.set_title(f'{pname}: ' + r'$\theta$ vs $|p|$', fontsize=10)

        # col 2: theta histogram
        ax = axes[row, 2]
        ax.hist(theta, bins=60, range=(0, thmax),
                histtype='stepfilled', color='lightskyblue', edgecolor='navy')
        ax.set_xlabel(r'$\theta$ (deg)', fontsize=10)
        ax.set_ylabel('Events', fontsize=9)
        ax.set_title(f'{pname} polar angle', fontsize=10)
        annotate(ax, theta)

    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f'Saved: {outfile}')

# ── main ───────────────────────────────────────────────────────────────────

print(f"Reading: {args.lund}")
events = read_lund(args.lund, args.nev)
print(f"  Read {len(events)} events")

if len(events) == 0:
    print("No events found — check file format.")
    exit(1)

d = extract(events)

# ── apply cuts ─────────────────────────────────────────────────────────────
mask = np.ones(len(events), dtype=bool)
if args.q2min > 0: mask &= d['Q2'] >= args.q2min
if args.q2max > 0: mask &= d['Q2'] <= args.q2max
if args.wmin  > 0: mask &= d['W']  >= args.wmin
if mask.sum() < len(events):
    print(f"  Cuts: {mask.sum()} / {len(events)} events pass "
          f"(Q2>={args.q2min}, Q2<={args.q2max}, W>={args.wmin})")
    d = {k: v[mask] if isinstance(v, np.ndarray) and v.ndim==1
            else tuple(x[mask] for x in v) if isinstance(v, tuple)
            else v
         for k, v in d.items()}

label    = args.label if args.label else os.path.basename(args.lund)
prefix   = args.out   if args.out   else os.path.splitext(args.lund)[0]

make_plot(d, f'Kinematics | {len(events)} events | {label}',
          prefix + '_validation.png')

make_kinematics_plot(d, f'Lab kinematics | {len(events)} events | {label}',
                     prefix + '_kinematics.png')
