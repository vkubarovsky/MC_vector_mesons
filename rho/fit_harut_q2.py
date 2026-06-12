"""
fit_harut_q2.py
Fit the Q² distribution of Harut sigma_T (narrow xB cut, NO t-cut) with A/(Q²)^{n/2}
and compare with:
  - DIFFRAD sigma_T only (sigT_born_events.lund)
  - Analytical DIFFRAD flux × sigma_T curve  (theory_q2_fit result: n=4.43)

Cuts: xB=[0.24,0.26], y=[0.10,0.80], W>=2.0, Q²=[1.0,5.5]  (NO t-cut)
Fit range: Q²=[1.2, 3.6] GeV²
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
SIGT_LUND  = os.path.join(HERE, 'sigT_born_events.lund')
HARUT_LUND = '/Users/vpk/Downloads/Harut_rho_MC_generator/Harut_sigma_T.lund'

cuts = dict(q2min=1.0, q2max=5.5, xbmin=0.24, xbmax=0.26,
            ymin=0.10, ymax=0.80, wmin=2.0)

def apply_cuts(evs, c):
    Q2 = np.array([e['Q2'] for e in evs])
    xB = np.array([e['xB'] for e in evs])
    y  = np.array([e['y']  for e in evs])
    W  = np.array([e['W']  for e in evs])
    m  = ((Q2 >= c['q2min']) & (Q2 <= c['q2max']) &
          (xB >= c['xbmin']) & (xB <= c['xbmax']) &
          (y  >= c['ymin'])  & (y  <= c['ymax'])  &
          (W  >= c['wmin']))
    return Q2[m]

print("Reading DIFFRAD sigma_T only ...", flush=True)
Q2_S = apply_cuts(read_lund_events(SIGT_LUND), cuts)
print(f"  N = {len(Q2_S)}")

print("Reading Harut sigma_T ...", flush=True)
Q2_H = apply_cuts(read_lund_events(HARUT_LUND), cuts)
print(f"  N = {len(Q2_H)}")

# ── Analytical theory: flux × sigma_T from DIFFRAD ───────────────────────
Mp, E, xB0, t0, bt = 0.938272, 10.2, 0.25, -0.5, 2.75
s = 2.0 * Mp * E
alpha = 1.0 / 137.036

def eps_f(Q2):
    y  = Q2 / (s * xB0)
    g2 = 4.0 * Mp**2 * xB0**2 / Q2
    num = 1 - y - 0.25*y**2*g2
    den = 1 - y + 0.50*y**2 + 0.25*y**2*g2
    return np.where((num > 0) & (den > 0), num/den, np.nan)

def theory_curve(Q2):
    y   = Q2 / (s * xB0)
    ep  = eps_f(Q2)
    fl  = (alpha/(2*np.pi)) * y**2/(1-ep) * (1-xB0)/(xB0*Q2)
    sT  = 30.0 * np.exp(bt*t0) * xB0/(1-xB0) / Q2
    return fl * sT

q2_fine = np.linspace(1.001, 4.49, 2000)
y_fine  = q2_fine / (s * xB0)
ok      = (y_fine > 0.10) & (y_fine < 0.80)
f_th    = np.where(ok, theory_curve(q2_fine), np.nan)

# ── Histogram setup ───────────────────────────────────────────────────────
NBINS    = 15
FIT_LO, FIT_HI = 1.2, 3.6
edges    = np.linspace(1.0, 4.5, NBINS + 1)
centers  = 0.5 * (edges[:-1] + edges[1:])
dQ2      = edges[1] - edges[0]

cnt_S, _ = np.histogram(Q2_S, bins=edges)
cnt_H, _ = np.histogram(Q2_H, bins=edges)

# ── Power-law fit ─────────────────────────────────────────────────────────
def model(Q2, A, n):
    return A * Q2**(-n / 2.0)

def do_fit(centers, counts, label):
    mask = (centers >= FIT_LO) & (centers <= FIT_HI) & (counts > 0)
    q2c  = centers[mask]
    nc   = counts[mask].astype(float)
    errs = np.maximum(np.sqrt(nc), 1.0)
    p0   = [nc[0] * q2c[0]**1.5, 3.0]
    popt, pcov = curve_fit(model, q2c, nc, sigma=errs, p0=p0,
                           bounds=([0, 0.5], [1e8, 10.0]))
    perr = np.sqrt(np.diag(pcov))
    chi2 = (((nc - model(q2c, *popt))**2) / np.maximum(nc, 1)).sum()
    ndf  = mask.sum() - 2
    print(f"  {label:<30}  n = {popt[1]:.3f} ± {perr[1]:.3f}"
          f"   χ²/ndf = {chi2/ndf:.2f}")
    return popt, perr

# Fit the analytical theory curve in the same range
th_mask  = ok & (q2_fine >= FIT_LO) & (q2_fine <= FIT_HI) & np.isfinite(f_th)
popt_th, pcov_th = curve_fit(model, q2_fine[th_mask], f_th[th_mask],
                              p0=[f_th[th_mask][0]*q2_fine[th_mask][0]**1.5, 4.0])
perr_th = np.sqrt(np.diag(pcov_th))
n_theory = popt_th[1]
print(f"\nFit  A/(Q²)^{{n/2}}  in Q²=[{FIT_LO},{FIT_HI}] GeV²:")
print(f"  {'DIFFRAD theory (analytic)':<30}  n = {n_theory:.3f} ± {perr_th[1]:.3f}")
popt_S, perr_S = do_fit(centers, cnt_S, 'DIFFRAD σT only (MC)')
popt_H, perr_H = do_fit(centers, cnt_H, 'Harut σT (MC)')

delta_n = popt_H[1] - popt_S[1]
print(f"\n  Harut − DIFFRAD σT:  Δn = {delta_n:.3f}")
print(f"  Harut − Theory:      Δn = {popt_H[1] - n_theory:.3f}")
print(f"  DIFFRAD − Theory:    Δn = {popt_S[1] - n_theory:.3f}")

# ── Normalise theory curve to histogram area ──────────────────────────────
def norm_to_hist(q2_arr, f_arr, ok_mask, counts):
    from scipy.integrate import trapezoid
    area_th   = trapezoid(f_arr[ok_mask], q2_arr[ok_mask])
    area_hist = counts.sum() * dQ2
    return f_arr * (area_hist / area_th)

th_S = norm_to_hist(q2_fine, f_th, ok, cnt_S)
th_H = norm_to_hist(q2_fine, f_th, ok, cnt_H)

# ── Plot ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle(r'Q² fits: $A/(Q^2)^{n/2}$  |  $x_B$∈[0.24,0.26], '
             r'no $t$-cut, $y$∈[0.10,0.80], $W$≥2.0  |  E=10.2 GeV',
             fontsize=11)

fit_colors = {'DIFFRAD σT': 'tab:orange', 'Harut σT': 'tab:green'}

for ax, (Q2_data, cnt, popt, perr, label, col, th_curve) in zip(axes, [
    (Q2_S, cnt_S, popt_S, perr_S, 'DIFFRAD σT only', 'tab:orange', th_S),
    (Q2_H, cnt_H, popt_H, perr_H, 'Harut σT',        'tab:green',  th_H),
]):
    errs = np.maximum(np.sqrt(cnt), 1.0)

    # data
    ax.errorbar(centers, cnt, yerr=errs, fmt='o', color=col,
                markersize=5, label=f'{label}  (N={len(Q2_data)})')

    # DIFFRAD theory curve (normalised)
    ax.plot(q2_fine[ok], th_curve[ok], 'k-', lw=1.8, alpha=0.7,
            label=fr'DIFFRAD flux×$\sigma_T$ theory  ($n_{{th}}={n_theory:.2f}$)')

    # power-law fit
    q2_fit_line = np.linspace(FIT_LO, FIT_HI, 400)
    ax.plot(q2_fit_line, model(q2_fit_line, *popt), '--', color=col, lw=2.2,
            label=fr'Fit $A/(Q^2)^{{n/2}}$:  $n = {popt[1]:.3f} \pm {perr[1]:.3f}$')

    ax.axvspan(FIT_LO, FIT_HI, alpha=0.07, color='red', label=f'Fit range [{FIT_LO},{FIT_HI}]')
    ax.set_yscale('log')
    ax.set_xlabel(r'$Q^2$ (GeV²)', fontsize=12)
    ax.set_ylabel('Events / bin', fontsize=11)
    ax.set_title(label, fontsize=11)
    ax.set_xlim(1.0, 4.5)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

fig.tight_layout()
out = os.path.join(HERE, 'fit_harut_q2_noTcut.png')
fig.savefig(out, dpi=150)
print(f"\nSaved: {out}")
