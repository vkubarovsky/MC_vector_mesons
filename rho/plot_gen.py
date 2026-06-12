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
                       (cuts applied ONLY if this argument is given)
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

cuts_cfg = {}
if args.input:
    if os.path.exists(args.input):
        cuts_cfg = read_gen_input(args.input)
        print(f"Cuts from: {args.input}")
    else:
        print(f"WARNING: input file {args.input} not found — no cuts applied")
else:
    print("No -input given — no kinematic cuts applied")

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
                 mlabel=r'$M_{h^+h^-}$',      mtitle=r'$J/\psi$ invariant mass',
                 dnames=(r'$h^+$', r'$h^-$')),
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

    Q2    = np.array([ev['Q2'] for ev in events])
    y     = np.array([ev['y']  for ev in events])
    W     = np.array([ev['W']  for ev in events])
    nu    = np.array([ev['nu'] for ev in events])
    t     = np.array([ev['t']  for ev in events])
    xB    = np.array([ev['xB'] for ev in events])
    ebeam = np.array([ev['ebeam'] for ev in events])

    ph_4v = pip + pim
    Mmeson = np.sqrt(np.maximum(0, m2v(ph_4v)))
    pT2    = ph_4v[:,1]**2 + ph_4v[:,2]**2
    Mmiss  = np.sqrt(np.maximum(0, m2v(pp + kgam)))
    E_miss = nu + Mp - pip[:,0] - pim[:,0] - pp[:,0] - kgam[:,0]
    Egam   = kgam[has_photon, 0]

    return dict(Q2=Q2, y=y, t=t, W=W, Mmeson=Mmeson, Mmiss=Mmiss,
                E_miss=E_miss, Egam=Egam, has_photon=has_photon, N=len(events),
                xB=xB, pT2=pT2, nu=nu, ebeam=ebeam,
                k2_kin=pkin(k2), pp_kin=pkin(pp),
                pip_kin=pkin(pip), pim_kin=pkin(pim))

# ── Power-law fit for Q² panel ───────────────────────────────────────────────
def fit_q2_powerlaw(ax, edges, counts, fit_range=(1.2, 3.6), color='red'):
    """Fit dN/dQ² = A*(Q²)^{-n/2} in fit_range and overlay the result."""
    centers = 0.5 * (edges[:-1] + edges[1:])
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
                      f'\n  range [{fit_range[0]}, {fit_range[1]}] GeV²')
        ax.legend(fontsize=8, loc='upper right')
    except Exception:
        pass

# ── Breit-Wigner helpers (rho, identical formula in diffrad_akushevich/harut) ─
def bw_rho_dNdM(M, M0=0.77526, G0=0.1502, mpi=0.13957):
    """Running-width Breit-Wigner dN/dM (sample_bw shape before sigma_T reweighting)."""
    p0   = np.sqrt(max(0.0, M0**2/4.0 - mpi**2))
    M2   = M**2
    pm   = np.sqrt(np.maximum(0.0, M2/4.0 - mpi**2))
    Grun = G0 * (pm/p0)**3 * (M0/M)
    return 2.0*M * Grun / ((M2 - M0**2)**2 + M0**2*Grun**2)

def sigt_difflt(M, q2, w2, t, tslope=3.0):
    """σ_T from difflt (Akushevich VMD), vectorized over M array."""
    M    = np.atleast_1d(np.asarray(M, dtype=float))
    amp  = 0.938272; amp2 = amp**2; p02 = 0.5; al_s = 0.25
    alpha_em = 1.0/137.036; pi_val = np.pi
    ggam_rho = 6.77e-6

    asx  = w2 + q2 - amp2
    aanu = asx / (2.0*amp)
    asxt = asx + t
    aeh  = asxt / (2.0*amp)
    amv2 = M**2
    aff2 = np.exp(tslope * t)
    atqt = t + q2 - amv2
    apt2 = (-(4.0*(aanu**2 + q2)*amv2 + 4.0*aanu*aeh*atqt
              - 4.0*aeh**2*q2 + atqt**2)) / (4.0*(aanu**2 + q2))
    apt2 = np.maximum(apt2, 0.0)

    axsb = (q2 + amv2 + apt2) / w2
    aq2b = (q2 + amv2 + apt2) / 4.0

    afm = np.where(
        apt2 <= p02,
        np.log((4.0*aq2b - apt2 + p02) / (apt2 + p02)),
        np.log((apt2 + p02) / (4.0*aq2b - apt2 + p02)
               * 4.0*aq2b**2 / np.maximum(apt2**2, 1e-30))
    )
    axsbgm = 3.0*(1.0 - axsb)**5
    denom  = 2.0*aq2b*(2.0*aq2b - apt2)*np.log(8.0*aq2b/p02)
    denom  = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    ask    = axsbgm * afm / denom

    sigt = al_s**2 * ggam_rho * amv2*M / 3.0 / alpha_em * pi_val**3 * ask**2 * aff2
    sigt = np.where((w2 < (amp + M)**2) | (axsb >= 1.0), 0.0, sigt)
    return sigt

def _harut_kin_accept(M_val, q2, w2, t_arr, ebeam, xB):
    """Fraction of events that pass M-dependent kinematic checks from diffrad_harut."""
    Mp_ = 0.938272; Mp2_ = Mp_**2
    amv2 = M_val**2
    ok = (np.sqrt(w2) >= Mp_ + M_val)
    tt1 = w2 - q2 - Mp2_
    tt2 = w2 - Mp2_ + amv2
    disc1 = tt1**2 + 4*q2*w2
    disc2 = tt2**2 - 4*amv2*w2
    ok &= (disc1 >= 0) & (disc2 >= 0)
    safe = ok.copy()
    sd1 = np.sqrt(np.maximum(0, disc1))
    sd2 = np.sqrt(np.maximum(0, disc2))
    tdmink = np.where(safe, -q2 + amv2 - 0.5/w2*(tt1*tt2 + sd1*sd2), 0)
    tdmaxk = np.where(safe, -q2 + amv2 - 0.5/w2*(tt1*tt2 - sd1*sd2), 0)
    ok &= (t_arr >= tdmink) & (t_arr <= tdmaxk)
    return ok


def bw_rho_overlay(ax, Mdata, mrange, nbins, Q2=None, W=None, t=None,
                   version='', ebeam=None, xB=None):
    """Overlay rho BW curve normalized to the histogram event count.
    For Harut: BW × kinematic_acceptance(M) with per-event ammax/Z_j weighting.
    For Akushevich: BW × <σ_T>(M) from difflt."""
    M0_ = 0.77526; G0_ = 0.1502; mpi_ = 0.13957; Mp_ = 0.938272
    ammin_ = 2*mpi_
    Mmin, Mmax = mrange
    bw_width   = (Mmax - Mmin) / nbins
    M_line     = np.linspace(Mmin, Mmax, 500)
    bw_vals    = bw_rho_dNdM(M_line)
    integral   = np.trapezoid(bw_vals, M_line)
    if integral <= 0:
        return
    N_in  = np.sum((Mdata >= Mmin) & (Mdata <= Mmax))
    scale = N_in * bw_width / integral
    ax.plot(M_line, scale * bw_vals, '--', color='forestgreen', lw=1.5, alpha=0.7,
            label=r'BW (running $\Gamma$)')

    is_harut = 'harut' in version.lower() if version else False

    if is_harut and Q2 is not None and W is not None and t is not None and len(Q2) > 0:
        N_ev = len(Q2)
        W2 = W**2
        ammax = np.minimum(W - Mp_, M0_ + 5*G0_)
        Z_full = np.trapezoid(bw_rho_dNdM(np.linspace(ammin_, M0_+5*G0_, 500)),
                              np.linspace(ammin_, M0_+5*G0_, 500))
        unique_am = np.unique(np.round(ammax, 5))
        norm_lut = {}
        for am in unique_am:
            if am <= ammin_:
                norm_lut[am] = 0.0
            else:
                mg = np.linspace(ammin_, am, 500)
                norm_lut[am] = np.trapezoid(bw_rho_dNdM(mg), mg)
        Z_per = np.array([norm_lut[round(am, 5)] for am in ammax])
        Z_per = np.where(Z_per > 0, Z_per, 1e-30)

        eff_bw = np.zeros_like(M_line)
        for i, M_val in enumerate(M_line):
            can_produce = ammax > M_val
            acc = _harut_kin_accept(M_val, Q2, W2, t, ebeam, xB) & can_produce
            eff_bw[i] = bw_rho_dNdM(M_val) * (acc.astype(float) / Z_per).sum()

        integral_h = np.trapezoid(eff_bw, M_line)
        if integral_h > 0:
            scale_h = N_in * bw_width / integral_h
            ax.plot(M_line, scale_h * eff_bw, '-', color='forestgreen', lw=2.5,
                    label=r'BW$\times$kin. accept.')

    elif Q2 is not None and W is not None and t is not None and len(Q2) > 0:
        sigt_per_evt = sigt_difflt(Mdata, Q2, W**2, t)
        n_mbins = 80
        m_edges = np.linspace(Mmin, Mmax, n_mbins + 1)
        m_centers = 0.5 * (m_edges[:-1] + m_edges[1:])
        mean_sigt = np.zeros(n_mbins)
        for i in range(n_mbins):
            mask = (Mdata >= m_edges[i]) & (Mdata < m_edges[i+1])
            if mask.sum() > 0:
                mean_sigt[i] = np.mean(sigt_per_evt[mask])
        from scipy.interpolate import interp1d
        good = mean_sigt > 0
        if good.sum() >= 2:
            f_sigt = interp1d(m_centers[good], mean_sigt[good],
                              kind='linear', fill_value='extrapolate')
            sigt_smooth = np.maximum(f_sigt(M_line), 0.0)
            bw_st = bw_vals * sigt_smooth
            integral_st = np.trapezoid(bw_st, M_line)
            if integral_st > 0:
                scale_st = N_in * bw_width / integral_st
                ax.plot(M_line, scale_st * bw_st, '-', color='forestgreen', lw=2.5,
                        label=r'BW$\times\langle\sigma_T\rangle$(M)')

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
    fit_q2_powerlaw(ax, edges_q2, counts_q2, color='red')

    # (0,1) xB
    ax = axes[0,1]
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
    ax.hist(mt, bins=60, range=(tlo, thi), **bkw)
    ax.set_xlabel(r'$-t$ (GeV$^2$)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$t$ distribution (exponential)')
    ax.set_yscale('log')
    annotate(ax, mt, ' GeV²')

    # (0,3) y
    ax = axes[0,3]
    ax.hist(d['y'], bins=50, range=(0.0, 1.0), **bkw)
    ax.set_xlabel(r'$y$', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$y$ distribution')
    annotate(ax, d['y'])

    # (1,0) W
    ax = axes[1,0]
    ax.hist(d['W'], bins=50, **bkw)
    ax.set_xlabel(r'$W$ (GeV)', fontsize=11)
    ax.set_ylabel('Events'); ax.set_title(r'$W$ distribution')
    annotate(ax, d['W'], ' GeV')

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
    if particle == 'rho':
        bw_rho_overlay(ax, d['Mmeson'], cfg['mrange'], 60,
                           Q2=d['Q2'], W=d['W'], t=d['t'],
                           version=args.version,
                           ebeam=d['ebeam'], xB=d['xB'])
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
    annotate(ax, d['nu'], ' GeV')

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
if 'ymin'  in cuts_cfg: cuts &= d['y']  >= cuts_cfg['ymin']
if 'ymax'  in cuts_cfg: cuts &= d['y']  <= cuts_cfg['ymax']
if 'wmin'  in cuts_cfg: cuts &= d['W']  >= cuts_cfg['wmin']
if 'tmin'  in cuts_cfg and cuts_cfg['tmin'] > 0:
    cuts &= mt >= cuts_cfg['tmin']
if 'tmax'  in cuts_cfg: cuts &= mt      <= cuts_cfg['tmax']

if cuts.sum() < N:
    print(f"  Cuts: {cuts.sum()} / {N} events pass ({100*cuts.sum()/N:.1f}%)")
    for key in ('Q2','xB','y','W','nu','t','Mmeson','Mmiss','E_miss','pT2','ebeam'):
        d[key] = d[key][cuts]
    for key in ('k2_kin','pp_kin','pip_kin','pim_kin'):
        d[key] = tuple(arr[cuts] for arr in d[key])
    d['has_photon'] = d['has_photon'][cuts]
    d['Egam']       = d['Egam'][d['has_photon']] if d['has_photon'].sum() > 0 else np.array([])
    d['N']          = int(cuts.sum())
    N               = d['N']
    n_hard          = d['has_photon'].sum()

if N == 0:
    print("ERROR: 0 events pass cuts — nothing to plot. Check cuts or LUND file kinematics.")
    exit(1)

print(f"  E_missing: mean={d['E_miss'].mean():.4f}, std={d['E_miss'].std():.4f} GeV")

# Load vdist alongside LUND or from job subdirectories matching input filename
vdist_file  = _base + '_vdist.dat'
lund_base   = os.path.basename(lund_file)  # e.g. rc_events.lund
vdist_name  = lund_base[:-5] + '_vdist.dat' if lund_base.endswith('.lund') else lund_base + '_vdist.dat'
vdat_files  = sorted(glob.glob(f'{_lund_dir}/job_*/{vdist_name}'))
if vdat_files:
    parts = [np.loadtxt(f) for f in vdat_files if os.path.getsize(f) > 0]
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
