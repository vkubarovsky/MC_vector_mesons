#!/usr/bin/env python3
"""
plot_gen.py  --  Validation plots for a single diffrad_gen LUND file.

Two output pages (names do NOT conflict with plot_compare.py):
  *_gen_val.png  --  3x4 grid of 1-D kinematic histograms
  *_gen_kin.png  --  4x3 grid: |p|, theta vs |p|, theta for each particle

Usage:
    python plot_gen.py [-lund FILE] [-input gen_input_born.dat]
                       [-particle rho|phi|jpsi] [-nev N]

Options:
    -lund FILE         LUND file to read (default: lund/born_events.lund)
    -input FILE        Generator input file for kinematic cuts
                       (auto-detected: gen_input_born.dat or gen_input_rc.dat)
    -particle NAME     rho, phi, or jpsi  (auto-detected from cwd if omitted)
    -nev N             max events to read (0 = all, default)

Cuts applied from input file: q2min, q2max, ymin, ymax, xbmin, xbmax, wmin, tmax

Output (next to the LUND file, names do not conflict with plot_compare.py):
    *_gen_val.png      3x4 grid of 1-D kinematic histograms
    *_gen_kin.png      4x3 grid: |p|, theta vs |p|, theta for each particle
"""

import numpy as np
import matplotlib.pyplot as plt
import os, glob, argparse
from scipy.optimize import curve_fit
from lund_reader import read_lund_events, resolve_lund_dir

parser = argparse.ArgumentParser(description='Validation plots for a single LUND file')
parser.add_argument('-version',  default='', help='Model version (e.g. akushevich, harut) — sets lund_VERSION as default lund dir')
parser.add_argument('-lund',     default=None,
                    help='LUND file (default: lund_VERSION/born_events.lund or lund/born_events.lund)')
parser.add_argument('-input',    default=None,
                    help='Generator input file for kinematic cuts (e.g. gen_input_born.dat)')
parser.add_argument('-particle', default='',
                    help='Particle: rho, phi, jpsi  (auto-detected from cwd if omitted)')
parser.add_argument('-nev',      default=0, type=int, help='Max events to read (0 = all)')
args      = parser.parse_args()
if args.lund:
    lund_file = args.lund
    _lund_dir = os.path.dirname(os.path.abspath(lund_file))
else:
    _lund_dir = resolve_lund_dir(f'lund_{args.version}') if args.version else 'lund'
    if not os.path.isdir(_lund_dir):
        import sys
        avail = [d for d in os.listdir('.') if d.startswith('lund_') and not d.endswith('.py') and (os.path.isdir(d) or os.path.isfile(d))]
        sys.exit(f"Error: lund directory '{_lund_dir}' not found.\n"
                 f"  Use -version to specify one of: {avail}\n"
                 f"  Example: python {os.path.basename(__file__)} -version akushevich")
    lund_file = f'{_lund_dir}/born_events.lund'
max_ev    = args.nev

# ── Parse generator input file for kinematic cuts ─────────────────────────────
def read_gen_input(fname):
    """Read gen_input_*.dat and return dict of kinematic limits."""
    keys = ['bmom','tmom','lepton','ivec','cutv','nev','seed',
            'q2min','q2max','ymin','ymax','tmin','tmax','tslope',
            'iborn','wmin','xbmin','xbmax']
    vals = {}
    with open(fname) as f:
        for i, line in enumerate(f):
            line = line.split('!')[0].strip()
            if line and i < len(keys):
                try: vals[keys[i]] = float(line)
                except ValueError: pass
    return vals

_input_file = args.input
if _input_file is None:
    # auto-detect: prefer born input, fall back to rc
    for candidate in ('gen_input_born.dat', 'gen_input_rc.dat'):
        if os.path.exists(candidate):
            _input_file = candidate
            break

cuts_cfg = {}
if _input_file and os.path.exists(_input_file):
    cuts_cfg = read_gen_input(_input_file)
    print(f"Cuts from: {_input_file}")
else:
    print("No input file found — no kinematic cuts applied")

# ── Particle-specific settings ────────────────────────────────────────────────
# Auto-detect from working directory name
_cwd = os.path.basename(os.getcwd()).lower()
particle = args.particle.lower() or (
    'jpsi' if 'jpsi' in _cwd else
    'phi'  if 'phi'  in _cwd else
    'rho')

PARTICLE_CFG = {
    'rho':  dict(label=r'$\rho^0$',  mass=0.7683,  mrange=(0.60, 0.95),
                 mlabel=r'$M_{\pi^+\pi^-}$', mtitle=r'$\rho^0$ invariant mass',
                 dnames=(r'$\pi^+$', r'$\pi^-$')),
    'phi':  dict(label=r'$\phi$',    mass=1.0195,  mrange=(0.98, 1.06),
                 mlabel=r'$M_{K^+K^-}$',      mtitle=r'$\phi$ invariant mass',
                 dnames=(r'$K^+$', r'$K^-$')),
    'jpsi': dict(label=r'$J/\psi$',  mass=3.0969,  mrange=(2.80, 3.30),
                 mlabel=r'$M_{e^+e^-}$',      mtitle=r'$J/\psi$ invariant mass',
                 dnames=(r'$e^+$', r'$e^-$')),
}
cfg = PARTICLE_CFG.get(particle, PARTICLE_CFG['rho'])
print(f"Particle: {particle}  ({cfg['label']})")

# Output PNG names (distinct from plot_compare.py)
_base          = lund_file[:-5] if lund_file.endswith('.lund') else lund_file
out_val        = _base + '_gen_val.png'
out_kin        = _base + '_gen_kin.png'

Mp = 0.938272

# ── Helpers ───────────────────────────────────────────────────────────────────

def m2v(p):
    return p[:,0]**2 - p[:,1]**2 - p[:,2]**2 - p[:,3]**2

def pkin(p4):
    pmag  = np.sqrt(p4[:,1]**2 + p4[:,2]**2 + p4[:,3]**2)
    theta = np.degrees(np.arctan2(np.sqrt(p4[:,1]**2 + p4[:,2]**2), p4[:,3]))
    phi   = np.degrees(np.arctan2(p4[:,2], p4[:,1]))
    return pmag, theta, phi

def annotate(ax, data, unit='', fmt='.4f', loc='right'):
    data = data[np.isfinite(data)]
    n    = len(data)
    mean = np.mean(data)
    sig  = np.std(data)
    w    = max(len(f'{mean:{fmt}}'), len(f'{sig:{fmt}}'), len(str(n)))
    txt  = '\n'.join([f'N     = {n:>{w}}',
                      f'Mean = {mean:>{w}{fmt}}{unit}',
                      f'Sigma = {sig:>{w}{fmt}}{unit}'])
    x, ha = (0.97, 'right') if loc == 'right' else (0.03, 'left')
    ax.text(x, 0.95, txt, transform=ax.transAxes,
            ha=ha, va='top', fontsize=9, fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow',
                      ec='gray', alpha=0.9))

def extract(events):
    has_photon = np.array([ev['has_photon'] for ev in events])
    k2   = np.array([ev['k2']   for ev in events])
    pip  = np.array([ev['pip']  for ev in events])
    pim  = np.array([ev['pim']  for ev in events])
    pp   = np.array([ev['pp']   for ev in events])
    kgam = np.array([ev['kgam'] for ev in events])

    Q2 = np.array([ev['Q2'] for ev in events])
    y  = np.array([ev['y']  for ev in events])
    W  = np.array([ev['W']  for ev in events])
    nu = np.array([ev['nu'] for ev in events])
    t  = np.array([ev['t']  for ev in events])
    xB = np.array([ev['xB'] for ev in events])

    ph_4v = pip + pim
    Mmeson = np.sqrt(np.maximum(0, m2v(ph_4v)))
    pT2    = ph_4v[:,1]**2 + ph_4v[:,2]**2
    Mmiss  = np.sqrt(np.maximum(0, m2v(pp + kgam)))
    E_miss = nu + Mp - pip[:,0] - pim[:,0] - pp[:,0] - kgam[:,0]
    Egam   = kgam[has_photon, 0]

    return dict(Q2=Q2, y=y, t=t, W=W, Mmeson=Mmeson, Mmiss=Mmiss,
                E_miss=E_miss, Egam=Egam, has_photon=has_photon, N=len(events),
                xB=xB, pT2=pT2, nu=nu,
                k2_kin=pkin(k2), pp_kin=pkin(pp),
                pip_kin=pkin(pip), pim_kin=pkin(pim))

# ── Page 1: 3x4 validation histograms ─────────────────────────────────────────

def make_validation(d, vv, title, outfile, facecolor, edgecolor):
    N      = d['N']
    n_hard = d['has_photon'].sum()
    bkw    = dict(histtype='stepfilled', color=facecolor, edgecolor=edgecolor)

    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    # (0,0) Q2
    ax = axes[0,0]
    counts_q2, edges_q2, _ = ax.hist(d['Q2'], bins=50, **bkw)
    ax.set_xlabel(r'$Q^2$ (GeV$^2$)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$Q^2$ distribution')
    ax.set_yscale('log')
    annotate(ax, d['Q2'], ' GeV²')
    # Fit A/(Q²)^(n/2) from first nonzero bin
    cen_q2 = 0.5*(edges_q2[:-1] + edges_q2[1:])
    nz = np.where(counts_q2 > 0)[0]
    if len(nz) >= 4:
        fr = (cen_q2[nz[0]], min(cen_q2[nz[-1]], 0.1))
        m = (cen_q2 >= fr[0]) & (cen_q2 <= fr[1]) & (counts_q2 > 0)
        q2c, nc = cen_q2[m], counts_q2[m].astype(float)
        errs = np.maximum(np.sqrt(nc), 1.0)
        try:
            def _pw(q2, A, n): return A * q2**(-n/2.0)
            p0 = [nc[0]*q2c[0]**(3.0/2.0), 3.0]
            po, pc = curve_fit(_pw, q2c, nc, sigma=errs, p0=p0,
                               bounds=([0,0.2],[1e9,12.0]))
            pe = np.sqrt(np.diag(pc))
            q2f = np.linspace(fr[0], fr[1], 300)
            ax.plot(q2f, _pw(q2f, *po), 'r-', lw=2,
                    label=fr'Fit $A/(Q^2)^{{n/2}}$: $n={po[1]:.2f}\pm{pe[1]:.2f}$'
                          f'\n  range [{fr[0]:.3f}, {fr[1]:.3f}] GeV²')
            ax.legend(fontsize=8, loc='upper right',
                      bbox_to_anchor=(0.97,0.75), bbox_transform=ax.transAxes)
        except Exception:
            pass

    # (0,1) xB
    ax = axes[0,1]
    if particle == 'jpsi':
        ax.hist(d['xB'], bins=50, range=(0.0, 0.016), **bkw)
    else:
        ax.hist(d['xB'], bins=50, **bkw)
    ax.set_xlabel(r'$x_B$', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$x_B$ distribution')
    ax.set_yscale('log')
    annotate(ax, d['xB'])

    # (0,2) t
    ax = axes[0,2]
    mt = -d['t']
    tlo = np.floor(np.min(mt) * 20) / 20
    thi = np.ceil (np.max(mt) * 20) / 20
    counts_t, edges_t, _ = ax.hist(mt, bins=60, range=(tlo, thi), **bkw)
    ax.set_xlabel(r'$-t$ (GeV$^2$)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$t$ distribution (exponential)')
    ax.set_yscale('log')
    annotate(ax, mt, ' GeV²')
    # Fit A·exp(-b·|t|) in range [1, 4]
    cen_t = 0.5*(edges_t[:-1] + edges_t[1:])
    mt_fit = (cen_t >= 1.0) & (cen_t <= 4.0) & (counts_t > 0)
    if mt_fit.sum() >= 4:
        tc, nc_t = cen_t[mt_fit], counts_t[mt_fit].astype(float)
        errs_t = np.maximum(np.sqrt(nc_t), 1.0)
        try:
            def _exp(x, A, b): return A * np.exp(-b * x)
            po_t, pc_t = curve_fit(_exp, tc, nc_t, sigma=errs_t,
                                   p0=[nc_t[0]*np.exp(tc[0]), 1.0])
            pe_t = np.sqrt(np.diag(pc_t))
            tf = np.linspace(1.0, 4.0, 200)
            ax.plot(tf, _exp(tf, *po_t), 'r-', lw=2,
                    label=fr'Fit $Ae^{{-b|t|}}$: $b={po_t[1]:.2f}\pm{pe_t[1]:.2f}$'
                          f'\n  range [1.0, 4.0] GeV²')
            ax.legend(fontsize=8, loc='upper right',
                      bbox_to_anchor=(0.97, 0.75), bbox_transform=ax.transAxes)
            print(f"  t-slope fit: b = {po_t[1]:.3f} +/- {pe_t[1]:.3f} GeV^-2")
        except Exception:
            pass

    # (0,3) y
    ax = axes[0,3]
    if particle == 'jpsi':
        ax.hist(d['y'], bins=50, range=(0.7, 1.0), **bkw)
        ax.set_xlim(0.7, 1.0)
    else:
        ax.hist(d['y'], bins=50, range=(0.0, 1.0), **bkw)
    ax.set_xlabel(r'$y$', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$y$ distribution')
    annotate(ax, d['y'], loc='left' if particle == 'jpsi' else 'right')

    # (1,0) W
    ax = axes[1,0]
    ax.hist(d['W'], bins=50, **bkw)
    ax.set_xlabel(r'$W$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$W$ distribution')
    annotate(ax, d['W'], ' GeV', loc='left' if particle == 'jpsi' else 'right')

    # (1,1) missing mass
    ax = axes[1,1]
    ax.hist(d['Mmiss'], bins=60, range=(0.85, 1.20), **bkw)
    ax.axvline(Mp, color='red', lw=1.5, ls='--', label=f'$M_p$={Mp:.4f} GeV')
    ax.set_xlabel(r'$M_x$ (GeV)', fontsize=10)
    ax.set_ylabel('Events'); ax.set_title(r'Missing mass $M_x$')
    ax.set_yscale('log')
    ax.legend(fontsize=9)
    annotate(ax, d['Mmiss'], ' GeV')

    # (1,2) pT2 meson
    ax = axes[1,2]
    ax.hist(d['pT2'], bins=50, **bkw)
    ax.set_xlabel(r'$p_T^2$ (GeV$^2/c^2$)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(cfg['label'] + r'  $p_T^2$ (lab)')
    ax.set_yscale('log')
    annotate(ax, d['pT2'], ' GeV²')

    # (1,3) meson invariant mass
    ax = axes[1,3]
    ax.hist(d['Mmeson'], bins=60, range=cfg['mrange'], **bkw)
    ax.axvline(cfg['mass'], color='red', lw=1.5, ls='--',
               label=f"{cfg['label']}={cfg['mass']:.4f} GeV")
    ax.set_xlabel(cfg['mlabel'] + ' (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(cfg['mtitle'])
    ax.legend(fontsize=9)
    annotate(ax, d['Mmeson'], ' GeV')

    # (2,0) v distribution
    ax = axes[2,0]
    if vv is not None and len(vv) > 0:
        ax.hist(vv, bins=40, **bkw)
        ax.set_xlabel(r'$v = M_X^2 - M_p^2$ (GeV$^2$)', fontsize=11)
        ax.set_ylabel('Hard events'); ax.set_title(r'Inelasticity $v$ distribution')
        ax.set_yscale('log')
        annotate(ax, vv, ' GeV²')
    else:
        ax.text(0.5, 0.5, 'No hard events\n(Born only)',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('Inelasticity v')

    # (2,1) photon energy
    ax = axes[2,1]
    if n_hard > 0:
        ax.hist(d['Egam'], bins=50, **bkw)
        ax.set_xlabel(r'$E_\gamma$ (GeV)', fontsize=11)
        ax.set_ylabel('Hard events'); ax.set_title(r'Radiated photon energy')
        ax.set_yscale('log')
        annotate(ax, d['Egam'], ' GeV')
    else:
        ax.text(0.5, 0.5, 'No hard events\n(Born only)',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('Radiated photon energy')

    # (2,2) missing energy
    ax = axes[2,2]
    ax.hist(d['E_miss'], bins=80, range=(-0.07, 0.07), **bkw)
    ax.set_xlabel(r'$E_{\rm miss} = E_{\rm init} - \sum E_{\rm final}$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'Missing energy (should be $\delta(0)$)')
    ax.set_yscale('log')
    annotate(ax, d['E_miss'], ' GeV')

    # (2,3) nu energy transfer
    ax = axes[2,3]
    ax.hist(d['nu'], bins=50, **bkw)
    ax.set_xlabel(r'$\nu$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'Energy transfer $\nu$')
    annotate(ax, d['nu'], ' GeV', loc='left' if particle == 'jpsi' else 'right')

    for ax in axes.flat:
        ax.tick_params(direction='in', which='both')
        ax.grid(True, alpha=0.4, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Saved: {outfile}")

# ── Page 2: 4x3 kinematics (|p|, theta vs |p|, theta) ─────────────────────────

def make_kinematics(d, title, outfile, facecolor, edgecolor):
    bkw = dict(histtype='stepfilled', color=facecolor, edgecolor=edgecolor)
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    plabels    = (r"$e'$", r'$p$', cfg['dnames'][0], cfg['dnames'][1])
    theta_maxs = (20, 90, 40, 40)
    kin_keys   = ['k2_kin', 'pp_kin', 'pip_kin', 'pim_kin']

    for row, (key, pname, thmax) in enumerate(zip(kin_keys, plabels, theta_maxs)):
        pmag, theta, _ = d[key]
        pmax = np.percentile(pmag, 99)

        # col 0: momentum magnitude
        ax = axes[row, 0]
        ax.hist(pmag, bins=60, **bkw)
        ax.set_xlabel(r'$|p|$ (GeV/c)', fontsize=10)
        ax.set_ylabel('Events')
        ax.set_title(f'{pname}  momentum')
        annotate(ax, pmag, ' GeV/c', fmt='.3f', loc='left' if row == 0 else 'right')

        # col 1: theta vs momentum (2D)
        ax = axes[row, 1]
        h = ax.hist2d(pmag, theta, bins=60, cmap='plasma', cmin=1,
                      range=[[0, pmax], [0, thmax]])
        plt.colorbar(h[3], ax=ax, pad=0.02)
        ax.set_xlabel(r'$|p|$ (GeV/c)', fontsize=10)
        ax.set_ylabel(r'$\theta$ (deg)', fontsize=10)
        ax.set_title(f'{pname}  '+r'$\theta$ vs $|p|$')

        # col 2: theta histogram
        ax = axes[row, 2]
        ax.hist(theta, bins=60, range=(0, thmax), **bkw)
        ax.set_xlabel(r'$\theta$ (deg)', fontsize=10)
        ax.set_ylabel('Events')
        ax.set_title(f'{pname}  polar angle')
        annotate(ax, theta[(theta >= 0) & (theta <= thmax)], ' deg', fmt='.2f')

    for ax in axes.flat:
        ax.tick_params(direction='in', which='both')
        ax.grid(True, alpha=0.4, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Saved: {outfile}")

# ── Read ───────────────────────────────────────────────────────────────────────

if not os.path.exists(lund_file):
    print(f"ERROR: {lund_file} not found."); exit(1)

print(f"Reading {lund_file} ...")
events = read_lund_events(lund_file, max_ev=max_ev)
N = len(events)
if N == 0:
    print(f"ERROR: no events read from {lund_file}"); exit(1)

d      = extract(events)
n_hard = d['has_photon'].sum()
print(f"  {N} events,  Hard: {n_hard} ({100*n_hard/N:.1f}%),  Soft: {N-n_hard}")

# ── Apply kinematic cuts from input file ──────────────────────────────────────
cuts = np.ones(N, dtype=bool)
mt   = -d['t']
# Only apply cuts that the generator actually enforces in sample_born:
# q2min/q2max, xbmin/xbmax, tmax, wmin.
# ymin/ymax are read by the generator but NOT applied — do not cut on them.
if 'q2min' in cuts_cfg: cuts &= d['Q2'] >= cuts_cfg['q2min']
if 'q2max' in cuts_cfg: cuts &= d['Q2'] <= cuts_cfg['q2max']
if 'xbmin' in cuts_cfg: cuts &= d['xB'] >= cuts_cfg['xbmin']
if 'xbmax' in cuts_cfg: cuts &= d['xB'] <= cuts_cfg['xbmax']
if 'wmin'  in cuts_cfg: cuts &= d['W']  >= cuts_cfg['wmin']
if 'tmin'  in cuts_cfg and cuts_cfg['tmin'] > 0:
    cuts &= mt >= cuts_cfg['tmin']
if 'tmax'  in cuts_cfg: cuts &= mt      <= cuts_cfg['tmax']

if cuts.sum() < N:
    print(f"  Cuts: {cuts.sum()} / {N} events pass ({100*cuts.sum()/N:.1f}%)")
    for key in ('Q2','xB','y','W','nu','t','Mmeson','Mmiss','E_miss','pT2'):
        d[key] = d[key][cuts]
    for key in ('k2_kin','pp_kin','pip_kin','pim_kin'):
        d[key] = tuple(arr[cuts] for arr in d[key])
    d['has_photon'] = d['has_photon'][cuts]
    d['Egam']       = d['Egam'][d['has_photon']] if d['has_photon'].sum() > 0 else np.array([])
    d['N']          = int(cuts.sum())
    N               = d['N']
    n_hard          = d['has_photon'].sum()

print(f"  E_missing: mean={d['E_miss'].mean():.4f}, std={d['E_miss'].std():.4f} GeV")

# Load vdist alongside LUND or from job subdirectories matching input filename
vdist_file  = _base + '_vdist.dat'
lund_base   = os.path.basename(lund_file)  # e.g. rc_events.lund
vdist_name  = lund_base[:-5] + '_vdist.dat' if lund_base.endswith('.lund') else lund_base + '_vdist.dat'
vdat_files  = sorted(glob.glob(f'{_lund_dir}/job_*/{vdist_name}'))
if vdat_files:
    parts = [np.atleast_1d(np.loadtxt(f)) for f in vdat_files if os.path.getsize(f) > 0]
    vv = np.concatenate(parts) if parts else np.array([])
    print(f"  vdist: {len(vv)} entries from {len(vdat_files)} job dirs")
elif os.path.exists(vdist_file) and os.path.getsize(vdist_file) > 0:
    vv = np.loadtxt(vdist_file)
    print(f"  vdist: {len(vv)} entries from {vdist_file}")
else:
    vv = np.array([])

# ── Titles ─────────────────────────────────────────────────────────────────────
fname_short = os.path.basename(lund_file)
pname       = cfg['label']
title_val   = (f'{fname_short}  |  {pname}  |  {N} events  |  '
               f'Hard: {n_hard} ({100*n_hard/N:.1f}%)  Soft: {N-n_hard}')
title_kin   = f'{fname_short}  |  {pname}  |  {N} events  |  Lab kinematics'

fc, ec = ('lightskyblue', 'navy') if 'born' in lund_file.lower() else ('lightsalmon', 'darkred')

# ── Plot ───────────────────────────────────────────────────────────────────────
make_validation(d, vv, title_val, out_val, fc, ec)
make_kinematics(d, title_kin, out_kin, fc, ec)

os.system(f'open {out_val} {out_kin}')
