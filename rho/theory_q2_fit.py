"""
theory_q2_fit.py
Compute the theoretical Q² shape:
    f(Q²) = Photon_flux(Q², xB, E) × sigma_T(Q², xB, t)
using DIFFRAD formulas, then fit A/(Q²)^{n/2} in Q²=[1.2, 3.6].

DIFFRAD formulas (from bornin + difflt in diffrad_gen.f90):
  flux   = (alpha/2pi) * y²/(1-eps) * (1-xB)/(xB*Q²)
  sigma_T = 30 * exp(bt*t) * xB/(1-xB) / Q²          [bt=2.75 GeV^{-2}]
  eps     = (1-y - y²*gamma²/4) / (1-y + y²/2 + y²*gamma²/4)
  gamma²  = 4*Mp²*xB²/Q²
  y       = Q²/(s*xB),   s = 2*Mp*E

Parameters fixed at xB=0.25, t=-0.5 GeV², E=10.2 GeV (current gen_input_born.dat).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ── Physics constants and fixed parameters ────────────────────────────────
Mp   = 0.938272          # proton mass, GeV
E    = 10.2              # beam energy, GeV
xB   = 0.25             # fixed xB
t    = -0.5             # fixed t, GeV²
bt   = 2.75             # sigma_T t-slope, GeV^{-2}
alpha = 1.0 / 137.036

s = 2.0 * Mp * E        # s ≈ 19.14 GeV²

print(f"Parameters: E={E} GeV,  xB={xB},  t={t} GeV²")
print(f"            s={s:.4f} GeV²,  bt={bt}")

# ── DIFFRAD formulas ──────────────────────────────────────────────────────
def eps(Q2):
    y      = Q2 / (s * xB)
    gamma2 = 4.0 * Mp**2 * xB**2 / Q2
    num    = 1.0 - y - 0.25 * y**2 * gamma2
    den    = 1.0 - y + 0.50 * y**2 + 0.25 * y**2 * gamma2
    return np.where((num > 0) & (den > 0), num / den, np.nan)

def photon_flux(Q2):
    y   = Q2 / (s * xB)
    ep  = eps(Q2)
    ok  = np.isfinite(ep) & (ep < 1) & (y > 0)
    fl  = np.where(ok,
          (alpha / (2.0 * np.pi)) * y**2 / (1.0 - ep) * (1.0 - xB) / (xB * Q2),
          0.0)
    return fl

def sigma_T(Q2):
    return 30.0 * np.exp(bt * t) * xB / (1.0 - xB) / Q2

def theory(Q2):
    return photon_flux(Q2) * sigma_T(Q2)

# ── Evaluate on a fine grid ───────────────────────────────────────────────
q2 = np.linspace(1.001, 5.5, 2000)
y_vals = q2 / (s * xB)
valid  = (y_vals > 0.10) & (y_vals < 0.80)   # y cuts from gen_input_born.dat

f = theory(q2)
f = np.where(valid, f, np.nan)

# ── Power-law fit in [1.2, 3.6] ──────────────────────────────────────────
FIT_LO, FIT_HI = 1.2, 3.6

fit_mask = valid & (q2 >= FIT_LO) & (q2 <= FIT_HI) & np.isfinite(f) & (f > 0)
q2_fit   = q2[fit_mask]
f_fit    = f[fit_mask]

def model(Q2, A, n):
    return A * Q2**(-n / 2.0)

popt, pcov = curve_fit(model, q2_fit, f_fit,
                       p0=[f_fit[0] * q2_fit[0]**1.5, 3.0],
                       bounds=([0, 0.5], [1e6, 10.0]))
perr = np.sqrt(np.diag(pcov))
A_fit, n_fit = popt

print(f"\nPower-law fit  A/(Q²)^{{n/2}}  in Q²=[{FIT_LO},{FIT_HI}] GeV²:")
print(f"  n = {n_fit:.4f} ± {perr[1]:.4f}")
print(f"  A = {A_fit:.4e} ± {perr[0]:.4e}")

# Residuals at key Q² values
print(f"\n  Q²   theory      fit       ratio")
for q2v in [1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 3.6]:
    th = theory(np.array([q2v]))[0]
    ft = model(q2v, *popt)
    if th > 0:
        print(f"  {q2v:.1f}  {th:.4e}  {ft:.4e}  {ft/th:.4f}")

# ── Plot ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left: theory + fit
ax = axes[0]
ax.semilogy(q2[valid], f[valid], 'b-', lw=2.5,
            label=r'$\Gamma(Q^2)\times\sigma_T(Q^2)$  DIFFRAD')
q2_line = np.linspace(FIT_LO, FIT_HI, 500)
ax.semilogy(q2_line, model(q2_line, *popt), 'r-', lw=2.0,
            label=fr'Fit $A/(Q^2)^{{n/2}}$: $n={n_fit:.3f}\pm{perr[1]:.3f}$'
                  f'\n  range [{FIT_LO}, {FIT_HI}] GeV²')
ax.axvspan(FIT_LO, FIT_HI, alpha=0.08, color='red', label='Fit range')
ax.set_xlabel(r'$Q^2$ (GeV²)', fontsize=12)
ax.set_ylabel(r'$\Gamma \times \sigma_T$  (arb. units)', fontsize=11)
ax.set_title(fr'DIFFRAD flux$\times\sigma_T$ at $x_B={xB}$, $t={t}$ GeV²,'
             f'  $E={E}$ GeV', fontsize=10)
ax.legend(fontsize=10)
ax.set_xlim(1.0, 5.0)
ax.grid(True, alpha=0.35)

# Right: ratio  theory / fit  (shows how well A/(Q²)^{n/2} describes the shape)
ax = axes[1]
ratio = f[valid] / model(q2[valid], *popt)
ax.plot(q2[valid], ratio, 'b-', lw=2)
ax.axvspan(FIT_LO, FIT_HI, alpha=0.08, color='red', label='Fit range')
ax.axhline(1.0, color='k', lw=1.0, ls='--')
ax.set_xlabel(r'$Q^2$ (GeV²)', fontsize=12)
ax.set_ylabel(r'Theory / fit', fontsize=11)
ax.set_title(fr'Ratio: $\Gamma\times\sigma_T$ / [$A/(Q^2)^{{n/2}}$]   '
             fr'($n={n_fit:.3f}$)', fontsize=10)
ax.set_xlim(1.0, 5.0)
ax.set_ylim(0, 2.5)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.35)

fig.tight_layout()
out = os.path.join(HERE, 'theory_q2_fit.png')
fig.savefig(out, dpi=150)
print(f"\nSaved: {out}")
