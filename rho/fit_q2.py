"""
fit_q2.py — fit Q² distributions with a physically motivated model and compare
with the analytical photon-flux × sigma_T(1/Q²) prediction.

The Q² distribution at fixed (xB, t) is:

  dN/dQ² ∝ flux(Q²) × sigma(Q²)

  flux(Q²) = (alpha/2pi) * y²/(1-eps) * (1-xB)/(xB*Q²)    [virtual photon flux]

With sigma_T ∝ (Q²)^{-n_T}  and  y = Q²/(s*xB):

  dN/dQ² ∝ Q^{4-n_T-2} / (1-eps(Q²))  =  Q^{2-n_T} / (1-eps(Q²))

For DIFFRAD model: n_T=2  →  dN/dQ² ∝ 1/(1-eps)         [theory reference]
For DIFFRAD + sigma_L (sigma_L ~ const, no 1/Q²):
  dN/dQ² ∝ [1 + eps * (sigma_L/sigma_T)*Q²] / (1-eps)    [rises vs pure sigT]
For Harut (unknown):
  fit free n_T to extract effective Q² power of sigma_T

Fit model: f(Q²; A, n_T) = A * Q^{2-n_T} / (1-eps(Q²))
  n_T=2 → pure sigma_T(1/Q²) from DIFFRAD model → reference
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lund_reader import read_lund_events

HERE = os.path.dirname(os.path.abspath(__file__))

DIFFRAD_LUND = os.path.expanduser('~/Downloads/DIFFRAD_lund/rho/born_events.lund')
SIGT_LUND    = os.path.join(HERE, 'sigT_born_events.lund')
HARUT_LUND   = '/Users/vpk/Downloads/Harut_rho_MC_generator/Harut_sigma_T.lund'

cuts = dict(q2min=1.0, q2max=5.5, xbmin=0.24, xbmax=0.26,
            ymin=0.10, ymax=0.80, wmin=2.0, tmin=0.4, tmax=0.6)

def apply_cuts(evs, c):
    Q2 = np.array([e['Q2'] for e in evs])
    xB = np.array([e['xB'] for e in evs])
    y  = np.array([e['y']  for e in evs])
    W  = np.array([e['W']  for e in evs])
    t  = np.array([e['t']  for e in evs])
    mt = -t
    mask = ((Q2 >= c['q2min']) & (Q2 <= c['q2max']) &
            (xB >= c['xbmin']) & (xB <= c['xbmax']) &
            (y  >= c['ymin'])  & (y  <= c['ymax'])  &
            (W  >= c['wmin'])  &
            (mt >= c['tmin'])  & (mt <= c['tmax']))
    return Q2[mask]

print("Reading DIFFRAD (sigT+sigL) ...", flush=True)
D = apply_cuts(read_lund_events(DIFFRAD_LUND), cuts)
print(f"  N={len(D)}")
print("Reading DIFFRAD (sigT only) ...", flush=True)
S = apply_cuts(read_lund_events(SIGT_LUND), cuts)
print(f"  N={len(S)}")
print("Reading Harut sigma_T ...", flush=True)
H = apply_cuts(read_lund_events(HARUT_LUND), cuts)
print(f"  N={len(H)}")

# ── Kinematics ────────────────────────────────────────────────────────────
Mp, E  = 0.938272, 10.2
s      = 2 * Mp * E     # ~19.14 GeV²
xB0    = 0.25           # representative xB
t0     = -0.5           # representative t (GeV²)
Q0     = 1.0            # reference scale for fit

def eps_exact(Q2, xB=xB0):
    y      = np.atleast_1d(Q2) / (s * xB)
    g2     = 4 * Mp**2 * xB**2 / np.atleast_1d(Q2)
    num    = 1 - y - 0.25 * y**2 * g2
    den    = 1 - y + 0.5  * y**2 + 0.25 * y**2 * g2
    return np.where((num > 0) & (den > 0), num / den, 0.0)

def y_ok(Q2, ymin=0.10, ymax=0.80, xB=xB0):
    y = np.atleast_1d(Q2) / (s * xB)
    return (y > ymin) & (y < ymax)

# ── Histogram setup ────────────────────────────────────────────────────────
NBINS     = 15
q2_edges  = np.linspace(1.0, 4.5, NBINS + 1)   # stop at y=0.8 limit
q2_cents  = 0.5 * (q2_edges[:-1] + q2_edges[1:])
dQ2       = q2_edges[1] - q2_edges[0]

# ── Analytical theory shapes ─────────────────────────────────────────────
q2_fine = np.linspace(1.001, 4.49, 2000)
eps_f   = eps_exact(q2_fine)
ok_f    = y_ok(q2_fine)

# shape for sigma_T ~ (Q²)^{-n_T}:  dN/dQ² ∝ Q^{2-n_T}/(1-eps)
def theory_shape(Q2, n_T):
    eps = eps_exact(Q2)
    ok  = y_ok(Q2)
    val = np.where(ok & (eps < 1), Q2**(2 - n_T) / (1 - eps), 0.0)
    return val

# DIFFRAD model: n_T=2
ref_shape = theory_shape(q2_fine, n_T=2.0)

# sigma_T + eps*sigma_L shape (DIFFRAD born):
bt, bl = 2.75, 4.25
sigT_c = 30 * np.exp(bt * t0) * xB0 / (1 - xB0)   # coeff of 1/Q²
sigL_c = 30 * np.exp(bl * t0) * xB0                 # const in Q²
R_L    = sigL_c / sigT_c    # sigma_L / (sigma_T * Q²) at Q2=1
# dN/dQ² ∝ [sigT(Q²) + eps*sigL] / Q² × y² / (1-eps)
#         = [sigT_c/Q² + eps*sigL_c] / Q² × Q⁴ / (1-eps)
#         = [sigT_c + eps*sigL_c*Q²] / (1-eps)
full_shape = np.where(ok_f & (eps_f < 1),
    (sigT_c + eps_f * sigL_c * q2_fine) / (1.0 - eps_f),
    0.0)

# ── Fit model: A × Q^{2-n_T} / (1-eps)  evaluated at bin centres ─────────
def fit_model(q2_arr, A, n_T):
    eps = eps_exact(q2_arr)
    ok  = y_ok(q2_arr)
    return np.where(ok & (eps < 1), A * q2_arr**(2 - n_T) / (1 - eps), 0.0)

def fit_full_model(q2_arr, A, n_T, r_L):
    """sigma_T(1/Q^n_T) + eps*sigma_L(1/Q^n_L) with r_L = sigma_L / sigma_T at Q2=1"""
    eps = eps_exact(q2_arr)
    ok  = y_ok(q2_arr)
    # dN/dQ² ∝ (sigT + eps*sigL) × y²/(1-eps)/Q²
    #         = [Q^{-n_T} + eps*r_L*Q^{-n_T}*1] × Q^{4-2}/(1-eps)
    # using our normalisation Q0=1: sigT ∝ Q^{-n_T}, sigL ∝ Q^{-n_T+n_diff}
    # Simplify: assume sigma_L/sigma_T = r_L × Q² (DIFFRAD-like: sigma_L~const, sigT~1/Q²)
    sigT = q2_arr**(-n_T)
    sigL = r_L * q2_arr**(1 - n_T)  # sigL ~ const when n_T=2: 1-2=-1 → Q^{-1}?
    # Actually for DIFFRAD: sigT~1/Q², sigL~const → ratio = sigL/sigT ~ Q²
    # Let's just use: dN/dQ² ∝ [1 + eps*r_L*Q²] * Q^{2-n_T}/(1-eps)
    val = np.where(ok & (eps < 1),
        A * (1.0 + eps * r_L * q2_arr) * q2_arr**(2 - n_T) / (1 - eps),
        0.0)
    return val

def do_fit(data, label, model='pure'):
    counts, _ = np.histogram(data, bins=q2_edges)
    errs      = np.maximum(np.sqrt(counts), 1.0)
    mask      = counts > 0
    if mask.sum() < 3:
        print(f"  {label}: too few bins to fit")
        return None, None, counts

    try:
        if model == 'pure':
            popt, pcov = curve_fit(fit_model, q2_cents[mask], counts[mask],
                                   sigma=errs[mask], p0=[counts.max(), 2.0],
                                   bounds=([0, 0.5], [1e7, 6.0]))
        else:  # full: fit A, n_T, r_L
            popt, pcov = curve_fit(fit_full_model, q2_cents[mask], counts[mask],
                                   sigma=errs[mask], p0=[counts.max(), 2.0, R_L],
                                   bounds=([0, 0.5, 0], [1e7, 6.0, 100]))
        perr = np.sqrt(np.diag(pcov))
        return popt, perr, counts
    except Exception as exc:
        print(f"  {label}: fit failed — {exc}")
        return None, None, counts

print("\nFits:  dN/dQ² = A × (Q²)^{2-n_T} / (1-ε)   [pure sigma_T model]")
print(f"  {'Sample':<30}  n_T   ±err    χ²/ndf   comment")
print(f"  {'-'*70}")

popt_D, perr_D, cnt_D = do_fit(D, 'DIFFRAD σT+εσL',  model='full')
popt_S, perr_S, cnt_S = do_fit(S, 'DIFFRAD σT only', model='pure')
popt_H, perr_H, cnt_H = do_fit(H, 'Harut σT',        model='pure')

# chi² / ndf helper
def chi2ndf(data, counts, popt, model='pure'):
    mask = counts > 0
    if popt is None: return float('nan'), 0
    pred = (fit_model(q2_cents, *popt) if model=='pure' else
            fit_full_model(q2_cents, *popt))
    c2 = (((counts[mask] - pred[mask])**2) / np.maximum(counts[mask], 1))[..., :]
    return c2.sum(), mask.sum() - len(popt)

# report
for label, popt, perr, counts, model in [
    ('DIFFRAD σT+εσL',  popt_D, perr_D, cnt_D, 'full'),
    ('DIFFRAD σT only', popt_S, perr_S, cnt_S, 'pure'),
    ('Harut σT',        popt_H, perr_H, cnt_H, 'pure'),
]:
    if popt is None: continue
    c2, ndf = chi2ndf(None, counts, popt, model)
    if model == 'full':
        print(f"  {label:<30}  n_T={popt[1]:.2f}±{perr[1]:.2f}  "
              f"r_L={popt[2]:.2f}±{perr[2]:.2f}  χ²/ndf={c2/max(ndf,1):.2f}")
    else:
        print(f"  {label:<30}  n_T={popt[1]:.2f}±{perr[1]:.2f}  "
              f"χ²/ndf={c2/max(ndf,1):.2f}")

print(f"\n  DIFFRAD model (reference):             n_T=2.00 (sigma_T ~ 1/Q²)")

# ── Plot ─────────────────────────────────────────────────────────────────
def normalise(shape, counts):
    """Scale shape to match total histogram area."""
    area_h = counts.sum() * dQ2
    # integrate over fine grid
    area_s = np.trapz(shape[ok_f], q2_fine[ok_f])
    return shape * (area_h / area_s) if area_s > 0 else shape

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
fig.suptitle(r'Q² fits: flux $\times$ $\sigma$ model  |  E=10.2 GeV  '
             r'  $x_B$∈[0.24,0.26]  $-t$∈[0.40,0.60]  y∈[0.10,0.80]  W≥2.0',
             fontsize=10)

panels = [
    (D, cnt_D, popt_D, perr_D, 'full',
     'DIFFRAD σT+εσL', 'tab:blue',
     normalise(full_shape, cnt_D), 'flux×(σT+εσL) ref (n_T=2)'),
    (S, cnt_S, popt_S, perr_S, 'pure',
     'DIFFRAD σT only', 'tab:orange',
     normalise(ref_shape, cnt_S), 'flux×σT ref (n_T=2)'),
    (H, cnt_H, popt_H, perr_H, 'pure',
     'Harut σT', 'tab:green',
     normalise(ref_shape, cnt_H), 'flux×σT ref (n_T=2)'),
]

for ax, (data, counts, popt, perr, model, label, col, th_curve, th_label) in \
        zip(axes, panels):

    errs = np.maximum(np.sqrt(counts), 1.0)

    # data histogram
    ax.errorbar(q2_cents, counts, yerr=errs, fmt='o', color=col,
                markersize=5, label=f'Data  (N={len(data)})')

    # theory reference (n_T=2 from DIFFRAD model)
    ax.plot(q2_fine[ok_f], th_curve[ok_f], 'k-', lw=2, alpha=0.7,
            label=th_label)

    # best fit
    if popt is not None:
        q2_plot = q2_fine[ok_f]
        if model == 'pure':
            fit_curve = fit_model(q2_plot, *popt)
            fit_label = f'Fit: n_T={popt[1]:.2f}±{perr[1]:.2f}'
        else:
            fit_curve = fit_full_model(q2_plot, *popt)
            fit_label = (f'Fit: n_T={popt[1]:.2f}±{perr[1]:.2f}, '
                         f'r_L={popt[2]:.2f}±{perr[2]:.2f}')
        ax.plot(q2_plot, fit_curve, '--', color=col, lw=2, label=fit_label)

    ax.set_xlabel(r'$Q^2$ (GeV²)', fontsize=12)
    ax.set_ylabel('Events / bin', fontsize=11)
    ax.set_title(label, fontsize=11)
    ax.set_xlim(1.0, 4.5)
    ax.legend(fontsize=8.5)
    ax.grid(True, alpha=0.35)

fig.tight_layout()
out = os.path.join(HERE, 'fit_q2.png')
fig.savefig(out, dpi=150)
print(f"\nSaved: {out}")

# ── Summary table ─────────────────────────────────────────────────────────
print(f"\n{'='*62}")
print(f"  n_T from  dN/dQ² = A × (Q²)^{{2-n_T}} / (1-ε):")
print(f"  {'Sample':<30}  n_T     expected")
rows = [
    ('DIFFRAD σT+εσL',  popt_D, 'n_eff < 2  (sigma_L rises Q²)'),
    ('DIFFRAD σT only', popt_S, '2.00  (sigma_T ~ 1/Q²)'),
    ('Harut σT',        popt_H, '2.00  (if same model)'),
]
for lbl, p, exp in rows:
    val = f'{p[1]:.2f}±{p[3]:.2f}' if (p is not None and len(p)>3) else \
          (f'{p[1]:.2f}' if p is not None else 'failed')
    if p is not None and len(p) >= 2:
        print(f"  {lbl:<30}  {p[1]:.2f}    {exp}")
    else:
        print(f"  {lbl:<30}  ---     {exp}")
