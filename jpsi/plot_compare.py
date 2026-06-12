#!/usr/bin/env python3
"""
plot_compare.py  --  Separate validation plots for Born and RC samples
Usage:  python plot_compare.py [-born born_events.lund] [-rc rc_events.lund]
"""

import numpy as np
import matplotlib.pyplot as plt
import os, glob, argparse
from scipy.optimize import curve_fit
from lund_reader import read_lund_events, resolve_lund_dir

parser = argparse.ArgumentParser(description='Validation plots for Born and RC samples')
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
                 f"  Example: python {os.path.basename(__file__)} -version akushevich")
    born_lund = args.born or f'{_lund_dir}/born_events.lund'
    rc_lund   = args.rc   or f'{_lund_dir}/rc_events.lund'
max_ev    = args.nev

def _png_name(lund):
    return lund[:-5]+'_validation.png' if lund.endswith('.lund') else lund+'_validation.png'
born_png = _png_name(born_lund)
rc_png   = _png_name(rc_lund)

def _kin_name(lund):
    return lund[:-5]+'_kinematics.png' if lund.endswith('.lund') else lund+'_kinematics.png'
born_kin = _kin_name(born_lund)
rc_kin   = _kin_name(rc_lund)

Mp    = 0.938272
Mjpsi = 3.09690

def read_lund(fname, max_ev=0):
    return read_lund_events(fname, max_ev=max_ev)

def m2v(p):
    """Vectorized m2 for array of shape (N,4)"""
    return p[:,0]**2 - p[:,1]**2 - p[:,2]**2 - p[:,3]**2

def pkin(p4):
    """Return (|p|, theta_deg, phi_deg) for array of 4-vectors (N,4)."""
    pmag  = np.sqrt(p4[:,1]**2 + p4[:,2]**2 + p4[:,3]**2)
    theta = np.degrees(np.arctan2(np.sqrt(p4[:,1]**2 + p4[:,2]**2), p4[:,3]))
    phi   = np.degrees(np.arctan2(p4[:,2], p4[:,1]))
    return pmag, theta, phi

def extract(events):
    N          = len(events)
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

    ph_4v     = pip + pim
    Mphi_reco = np.sqrt(np.maximum(0, m2v(ph_4v)))
    pT2       = ph_4v[:,1]**2 + ph_4v[:,2]**2
    Mmiss     = np.sqrt(np.maximum(0, m2v(pp + kgam)))
    E_miss    = nu + Mp - pip[:,0] - pim[:,0] - pp[:,0] - kgam[:,0]
    Egam      = kgam[has_photon, 0]

    return dict(Q2=Q2, y=y, t=t, W=W, Mphi=Mphi_reco, Mmiss=Mmiss,
                E_miss=E_miss, Egam=Egam, has_photon=has_photon, N=N,
                xB=xB, pT2=pT2, nu=nu,
                k2_kin=pkin(k2), pp_kin=pkin(pp),
                pip_kin=pkin(pip), pim_kin=pkin(pim))

# ── Stats annotation ──────────────────────────────────────────────────────────
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

def fit_q2_powerlaw(ax, edges, counts, fit_range=None, color='red'):
    """Fit A/(Q²)^(n/2) to the Q² histogram.
    If fit_range is None, auto-detect: first nonzero bin to last bin."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    if fit_range is None:
        nonzero = np.where(counts > 0)[0]
        if len(nonzero) < 4:
            return
        fit_range = (centers[nonzero[0]], min(centers[nonzero[-1]], 0.1))
    mask    = (centers >= fit_range[0]) & (centers <= fit_range[1]) & (counts > 0)
    if mask.sum() < 4:
        return
    q2c, nc = centers[mask], counts[mask].astype(float)
    errs    = np.maximum(np.sqrt(nc), 1.0)
    def model(q2, A, n):
        return A * q2**(-n / 2.0)
    try:
        p0   = [nc[0] * q2c[0]**(3.0 / 2.0), 3.0]
        popt, pcov = curve_fit(model, q2c, nc, sigma=errs, p0=p0,
                               bounds=([0, 0.2], [1e9, 12.0]))
        perr = np.sqrt(np.diag(pcov))
        A, n = popt
        q2_line = np.linspace(fit_range[0], fit_range[1], 300)
        ax.plot(q2_line, model(q2_line, A, n), '-', color=color, lw=2.0,
                label=fr'Fit $A/(Q^2)^{{n/2}}$: $n={n:.2f}\pm{perr[1]:.2f}$'
                      f'\n  range [{fit_range[0]:.3f}, {fit_range[1]:.3f}] GeV²')
        ax.legend(fontsize=8, loc='upper right',
                  bbox_to_anchor=(0.97, 0.75), bbox_transform=ax.transAxes)
    except Exception:
        pass

# ── Main plot function ────────────────────────────────────────────────────────
def make_plot(d, vv, title, outfile, facecolor, edgecolor):
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
    fit_q2_powerlaw(ax, edges_q2, counts_q2, color='red')

    # (0,1) xB
    ax = axes[0,1]
    ax.hist(d['xB'], bins=50, range=(0.0, 0.016), **bkw)
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
    ax.hist(d['y'], bins=50, range=(0.7, 1.0), **bkw)
    ax.set_xlim(0.7, 1.0)
    ax.set_xlabel(r'$y$', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$y$ distribution')
    annotate(ax, d['y'], loc='left')

    # (1,0) W
    ax = axes[1,0]
    ax.hist(d['W'], bins=50, **bkw)
    ax.set_xlabel(r'$W$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$W$ distribution')
    annotate(ax, d['W'], ' GeV', loc='left')

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
    ax.set_ylabel('Events'); ax.set_title(r'Meson $p_T^2$ (lab)')
    ax.set_yscale('log')
    annotate(ax, d['pT2'], ' GeV²')

    # (1,3) J/psi mass
    ax = axes[1,3]
    ax.hist(d['Mphi'], bins=60, range=(2.90, 3.30), **bkw)
    ax.axvline(Mjpsi, color='red', lw=1.5, ls='--', label=f'$M_{{J/\\psi}}$={Mjpsi:.4f}')
    ax.set_xlabel(r'$M_{e^+e^-}$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$J/\psi$ invariant mass')
    ax.legend(fontsize=9)
    annotate(ax, d['Mphi'], ' GeV')

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
    annotate(ax, d['nu'], ' GeV', loc='left')

    for ax in axes.flat:
        ax.tick_params(direction='in', which='both')
        ax.grid(True, alpha=0.4, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Saved: {outfile}")
    return fig

# ── Kinematics page: p, theta vs p, phi for each final-state particle ─────────
def make_kinematics_plot(d, title, outfile, facecolor, edgecolor, plabels,
                          theta_maxs=(20, 90, 40, 40)):
    bkw = dict(histtype='stepfilled', color=facecolor, edgecolor=edgecolor)
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    kin_keys = ['k2_kin', 'pp_kin', 'pip_kin', 'pim_kin']
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

        # col 1: theta vs momentum (2D histogram)
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
    return fig

# ── Read data ─────────────────────────────────────────────────────────────────
print(f"Reading {born_lund} ...")
born = extract(read_lund(born_lund, max_ev))
print(f"  {born['N']} events,  Hard: {born['has_photon'].sum()}")

print(f"Reading {rc_lund} ...")
rc = extract(read_lund(rc_lund, max_ev))
n_hard = rc['has_photon'].sum()
print(f"  {rc['N']} events,  Hard: {n_hard} ({100*n_hard/rc['N']:.1f}%)")

# Collect vdist from job subdirectories or alongside the rc lund file
vdist_file = rc_lund[:-5]+'_vdist.dat' if rc_lund.endswith('.lund') else rc_lund+'_vdist.dat'
vdat_files = sorted(glob.glob(f'{_lund_dir}/job_*/rc_events_vdist.dat'))
if vdat_files:
    parts = [np.loadtxt(f, ndmin=1) for f in vdat_files if os.path.getsize(f) > 0]
    vv = np.concatenate(parts) if parts else np.array([])
    print(f"  vdist: {len(vv)} entries from {len(vdat_files)} job dirs")
elif os.path.exists(vdist_file) and os.path.getsize(vdist_file) > 0:
    vv = np.loadtxt(vdist_file)
else:
    vv = np.array([])

# ── Generate plots ────────────────────────────────────────────────────────────
make_plot(born, None,
          f'Born MC  |  {born["N"]} events  |  No radiative corrections',
          born_png, 'lightskyblue', 'navy')

make_plot(rc, vv,
          f'RC MC  |  {rc["N"]} events  |  '
          f'Hard: {n_hard} ({100*n_hard/rc["N"]:.1f}%)  '
          f'Soft: {rc["N"]-n_hard}',
          rc_png, 'lightsalmon', 'darkred')

PLABELS = [r"$e'$", r'$p$', r'$e^+$', r'$e^-$']

make_kinematics_plot(born,
    f'Born MC  |  {born["N"]} events  |  Lab kinematics',
    born_kin, 'lightskyblue', 'navy', PLABELS)

os.system(f'open {born_png} {rc_png} {born_kin}')
