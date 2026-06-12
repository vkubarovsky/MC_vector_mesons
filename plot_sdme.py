#!/usr/bin/env python3
"""
Plot SDME angular distribution test results.
Reads .dat files from test_sdme_plots and generates publication-quality plots.
Output: PNG files for the LaTeX report.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['axes.labelsize'] = 14
matplotlib.rcParams['axes.titlesize'] = 13

outdir = '/Users/vpk/OneDrive/Work/2026_DIFFRAD/report/figures'
import os
os.makedirs(outdir, exist_ok=True)

datadir = '/Users/vpk/OneDrive/Work/2026_DIFFRAD/MC_vector_mesons'

def load(fname):
    d = np.loadtxt(os.path.join(datadir, fname))
    return d[:,0], d[:,1]

# ── Analytic curves ──────────────────────────────────────
def W0_costh(ct, r04_00):
    """W^0 integrated over phi and Phi (only r04_00 nonzero)"""
    return (3.0/4.0) * (0.5*(1-r04_00)*(1-ct**2) + r04_00*ct**2)

def schc_costh(ct):
    """SCHC+NPE: r04_00=0 -> sin^2 theta"""
    return (3.0/8.0) * (1 - ct**2)  # (3/4)*(1/2)*sin^2

# ═══════════════════════════════════════════════════════════
# Figure 1: Three scenarios side by side — cos(theta)
# ═══════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# Panel (a): Isotropic
ct, n = load('sdme_iso_costh.dat')
dct = ct[1] - ct[0]
norm = np.sum(n) * dct
axes[0].bar(ct, n/norm, width=dct*0.9, color='steelblue', alpha=0.7, label='MC')
axes[0].axhline(0.5, color='red', ls='--', lw=2, label=r'$W=\frac{1}{2}$ (flat)')
axes[0].set_xlabel(r'$\cos\vartheta$')
axes[0].set_ylabel(r'$W(\cos\vartheta)$ (normalized)')
axes[0].set_title(r'(a) Isotropic: all $r^\alpha_{ij}=0$')
axes[0].legend(fontsize=10)
axes[0].set_ylim(0, 1.0)

# Panel (b): SCHC+NPE
ct, n = load('sdme_schc_costh.dat')
norm = np.sum(n) * dct
ct_th = np.linspace(-1, 1, 200)
axes[1].bar(ct, n/norm, width=dct*0.9, color='steelblue', alpha=0.7, label='MC')
axes[1].plot(ct_th, schc_costh(ct_th), 'r-', lw=2,
             label=r'$\frac{3}{8}\sin^2\vartheta$')
axes[1].set_xlabel(r'$\cos\vartheta$')
axes[1].set_title(r'(b) SCHC+NPE: $r^{04}_{00}=0$')
axes[1].legend(fontsize=10)
axes[1].set_ylim(0, 0.55)

# Panel (c): r04_00 = 0.5
ct, n = load('sdme_r04_costh.dat')
norm = np.sum(n) * dct
axes[2].bar(ct, n/norm, width=dct*0.9, color='steelblue', alpha=0.7, label='MC')
axes[2].plot(ct_th, W0_costh(ct_th, 0.5), 'r-', lw=2,
             label=r'$\frac{3}{4}[\frac{1}{4}\sin^2\vartheta + \frac{1}{2}\cos^2\vartheta]$')
axes[2].set_xlabel(r'$\cos\vartheta$')
axes[2].set_title(r'(c) $r^{04}_{00}=0.5$ (T+L mix)')
axes[2].legend(fontsize=10, loc='lower center')
axes[2].set_ylim(0, 1.0)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'sdme_costheta.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved sdme_costheta.png')

# ═══════════════════════════════════════════════════════════
# Figure 2: phi and Phi for SCHC+NPE
# ═══════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# phi distribution
ph, n = load('sdme_schc_phi.dat')
dph = ph[1] - ph[0]
norm = np.sum(n) * dph
ph_th = np.linspace(0, 360, 200)
# Analytic: integrated over costh and Phi with SCHC+NPE SDMEs
# phi dependence from r04_1-1=0 -> no cos2phi in W0
# but from r1_1-1=0.5 via eps*cos2Phi -> after Phi integration, cancels
# Actually, integrating W over Phi and costh:
# only W0 survives Phi integration. W0 has r04_1-1=0 for SCHC+NPE
# so phi should be flat after Phi integration
# But we sample (costh, phi, Phi) jointly, so phi is NOT marginalised
# The phi dist includes correlation with Phi
axes[0].bar(ph, n/norm, width=dph*0.9, color='coral', alpha=0.7, label='MC')
axes[0].axhline(1.0/360.0, color='gray', ls=':', lw=1)
axes[0].set_xlabel(r'$\varphi$ (deg)')
axes[0].set_ylabel('Normalized counts')
axes[0].set_title(r'(a) $\varphi$ distribution (SCHC+NPE)')
axes[0].legend(fontsize=10)

# Phi distribution
PH, n = load('sdme_schc_Phi.dat')
dPH = PH[1] - PH[0]
norm = np.sum(n) * dPH
PH_th = np.linspace(0, 360, 200)
# Analytic for SCHC+NPE integrated over costh and phi:
# From W1 term: -eps*cos2Phi * (r1_11*(2/3) + r1_00*(1/3))
# r1_11 = 0.5, r1_00 = 0, so contribution = -eps*cos2Phi*(1/3)
# From W0: 1/(4pi) after integration
# W_Phi ~ 1/(2pi) * [1 - eps*(2/3)*r1_11*cos2Phi] ...
# Actually easier to just show MC vs a fitted curve
axes[1].bar(PH, n/norm, width=dPH*0.9, color='coral', alpha=0.7, label='MC')
# Fit: a + b*cos(2*Phi)
PH_rad = PH * np.pi / 180.0
from scipy.optimize import curve_fit
def cos2phi_func(x, a, b):
    return a + b * np.cos(2 * x * np.pi / 180.0)
popt, _ = curve_fit(cos2phi_func, PH, n/norm)
axes[1].plot(PH_th, cos2phi_func(PH_th, *popt), 'r-', lw=2,
             label=r'Fit: $a + b\cos 2\Phi$')
axes[1].set_xlabel(r'$\Phi$ (deg)')
axes[1].set_title(r'(b) $\Phi$ distribution (SCHC+NPE)')
axes[1].legend(fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'sdme_phi_angles.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved sdme_phi_angles.png')

# ═══════════════════════════════════════════════════════════
# Figure 3: HERMES-like realistic scenario — all 3 angles
# ═══════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ct, n = load('sdme_hermes_costh.dat')
dct = ct[1] - ct[0]
norm = np.sum(n) * dct
axes[0].bar(ct, n/norm, width=dct*0.9, color='forestgreen', alpha=0.7, label='MC')
ct_th = np.linspace(-1, 1, 200)
axes[0].plot(ct_th, W0_costh(ct_th, 0.3), 'r-', lw=2,
             label=r'$W^0$ with $r^{04}_{00}=0.3$')
axes[0].set_xlabel(r'$\cos\vartheta$')
axes[0].set_ylabel(r'$W(\cos\vartheta)$')
axes[0].set_title(r'(a) $\cos\vartheta$')
axes[0].legend(fontsize=10)
axes[0].set_ylim(0, 0.8)

ph, n = load('sdme_hermes_phi.dat')
dph = ph[1] - ph[0]
norm = np.sum(n) * dph
axes[1].bar(ph, n/norm, width=dph*0.9, color='forestgreen', alpha=0.7, label='MC')
axes[1].set_xlabel(r'$\varphi$ (deg)')
axes[1].set_title(r'(b) $\varphi$')
axes[1].legend(fontsize=10)

PH, n = load('sdme_hermes_Phi.dat')
dPH = PH[1] - PH[0]
norm = np.sum(n) * dPH
axes[2].bar(PH, n/norm, width=dPH*0.9, color='forestgreen', alpha=0.7, label='MC')
popt2, _ = curve_fit(cos2phi_func, PH, n/norm)
axes[2].plot(PH_th, cos2phi_func(PH_th, *popt2), 'r-', lw=2,
             label=r'Fit: $a + b\cos 2\Phi$')
axes[2].set_xlabel(r'$\Phi$ (deg)')
axes[2].set_title(r'(c) $\Phi$')
axes[2].legend(fontsize=10)

fig.suptitle(r'Realistic electroproduction: $r^{04}_{00}=0.3$, $r^1_{11}=0.4$, $r^1_{1-1}=0.3$, $r^5_{00}=0.15$, $\varepsilon=0.85$',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'sdme_hermes_like.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved sdme_hermes_like.png')

# ═══════════════════════════════════════════════════════════
# Figure 4: Comparison of cos(theta) for different r04_00
# ═══════════════════════════════════════════════════════════
fig, ax = plt.subplots(1, 1, figsize=(7, 5))

ct_th = np.linspace(-1, 1, 200)
for r00, color, label in [(0.0, 'blue', r'$r^{04}_{00}=0$ (pure $\sigma_T$)'),
                            (0.3, 'green', r'$r^{04}_{00}=0.3$'),
                            (0.5, 'orange', r'$r^{04}_{00}=0.5$'),
                            (1.0, 'red', r'$r^{04}_{00}=1.0$ (pure $\sigma_L$)')]:
    ax.plot(ct_th, W0_costh(ct_th, r00), '-', color=color, lw=2.5, label=label)

ax.set_xlabel(r'$\cos\vartheta$', fontsize=14)
ax.set_ylabel(r'$W^0(\cos\vartheta)$', fontsize=14)
ax.set_title(r'Decay angular distribution $W^0$ vs $r^{04}_{00}$')
ax.legend(fontsize=11)
ax.set_ylim(0, 1.0)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'sdme_r00_comparison.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved sdme_r00_comparison.png')

print('\nAll plots saved to:', outdir)
