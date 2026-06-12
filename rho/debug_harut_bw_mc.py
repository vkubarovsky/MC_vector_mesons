#!/usr/bin/env python3
"""
debug_harut_bw_mc.py -- Monte Carlo replication of the Harut BW shape.

For each event in the LUND file, replicate exactly what diffrad_harut does:
  1. Use the event's (Q2, xB, t, W, ebeam) kinematics
  2. Sample M from running-width BW with per-event ammax
  3. Apply all M-dependent kinematic checks from the main loop
  4. Histogram the accepted M values

This gives the TRUE predicted M distribution without analytical approximations.
Also try a simple BW fit with free parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from lund_reader import read_lund_events

M0   = 0.77526
G0   = 0.1502
mpi  = 0.13957
Mp   = 0.938272
Mp2  = Mp**2
ammin = 2*mpi
p0 = np.sqrt(max(0, M0**2/4 - mpi**2))

fname = '/Users/vpk/Downloads/DIFFRAD_lund/harut/rho/born_events.lund'
print(f"Reading {fname} ...")
events = read_lund_events(fname, max_ev=100000)
N = len(events)
print(f"  {N} events loaded")

Q2 = np.array([e['Q2'] for e in events])
W  = np.array([e['W']  for e in events])
W2 = W**2
t  = np.array([e['t']  for e in events])
xB = np.array([e['xB'] for e in events])
ebeam = np.array([e['ebeam'] for e in events])
Mmeson = np.array([np.sqrt(max(0, (e['pip'][0]+e['pim'][0])**2
                                 -(e['pip'][1]+e['pim'][1])**2
                                 -(e['pip'][2]+e['pim'][2])**2
                                 -(e['pip'][3]+e['pim'][3])**2))
                   for e in events])

s   = 2*Mp*ebeam + Mp2
ys  = Q2 / (s * xB)
x   = s * (1 - ys)
sx  = s - x
sxp = s + x
aly = sx**2 + 4*Mp2*Q2

ammax_per_evt = np.minimum(W - Mp, M0 + 5*G0)

# ── sample_bw in Python (vectorized for many samples per event) ──────────

def sample_bw_python(ammax_val, n_samples, rng):
    """Sample n_samples masses from running-width BW in [ammin, ammax_val]."""
    if ammax_val <= ammin:
        return np.array([])

    thmin = np.arctan((ammin**2 - M0**2) / (M0 * G0))
    thmax = np.arctan((ammax_val**2 - M0**2) / (M0 * G0))

    # Find ratio_max by scanning
    th_scan = np.linspace(thmin, thmax, 500)
    amv2_scan = M0**2 + M0*G0*np.tan(th_scan)
    amv_scan = np.sqrt(np.maximum(ammin**2, amv2_scan))
    pm_scan = np.sqrt(np.maximum(0, amv_scan**2/4 - mpi**2))
    gr_scan = G0*(pm_scan/p0)**3*(M0/amv_scan)
    bwr_scan = 1.0/((amv2_scan - M0**2)**2 + M0**2*gr_scan**2)
    bwc_scan = 1.0/((amv2_scan - M0**2)**2 + M0**2*G0**2)
    ratio_scan = bwr_scan/bwc_scan * (gr_scan/G0)
    ratio_max = max(1.0, ratio_scan.max())

    samples = []
    while len(samples) < n_samples:
        batch = max(n_samples * 3, 1000)
        theta = thmin + rng.random(batch) * (thmax - thmin)
        amv2 = M0**2 + M0*G0*np.tan(theta)
        amv = np.sqrt(np.maximum(ammin**2, amv2))
        pm = np.sqrt(np.maximum(0, amv**2/4 - mpi**2))
        gamrun = G0*(pm/p0)**3*(M0/amv)
        bw_run = 1.0/((amv2 - M0**2)**2 + M0**2*gamrun**2)
        bw_const = 1.0/((amv2 - M0**2)**2 + M0**2*G0**2)
        ratio = bw_run/bw_const * (gamrun/G0) / ratio_max
        accept = rng.random(batch) <= ratio
        samples.extend(amv[accept].tolist())
    return np.array(samples[:n_samples])


def check_all_fortran(M_arr, q2_val, w2_val, t_val, sx_val, sxp_val, aly_val, s_val, x_val):
    """Apply all M-dependent checks from diffrad_harut main loop."""
    amv2 = M_arr**2

    ok = (np.sqrt(w2_val) >= Mp + M_arr)

    tt1 = w2_val - q2_val - Mp2
    tt2 = w2_val - Mp2 + amv2
    disc1 = tt1**2 + 4*q2_val*w2_val
    disc2 = tt2**2 - 4*amv2*w2_val
    ok &= (disc1 >= 0) & (disc2 >= 0)

    sd1 = np.sqrt(np.maximum(0, disc1))
    sd2 = np.sqrt(np.maximum(0, disc2))
    tdmink = np.where(ok, -q2_val + amv2 - 0.5/w2_val*(tt1*tt2 + sd1*sd2), 0)
    tdmaxk = np.where(ok, -q2_val + amv2 - 0.5/w2_val*(tt1*tt2 - sd1*sd2), 0)
    ok &= (t_val >= tdmink) & (t_val <= tdmaxk)

    sxt = sx_val + t_val
    tq = q2_val + t_val - amv2
    vmax_kin = np.where(ok,
        tt2 + 0.5/q2_val*(-tt1*tq
            - np.sqrt(np.maximum(0, tt1**2 + 4*q2_val*w2_val))
            * np.sqrt(np.maximum(0, tq**2 + 4*amv2*q2_val)))
        - 1e-8, 1.0)
    ok &= (vmax_kin > 0)

    aa1 = np.where(ok,
        (q2_val*sxp_val*sxt - (s_val*sx_val + 2*Mp2*q2_val)*tq) / (2*aly_val), 0)
    aa2 = np.where(ok,
        (q2_val*sxp_val*sxt - (x_val*sx_val - 2*Mp2*q2_val)*tq) / (2*aly_val), 0)
    ssh = x_val + q2_val - aa2/2
    xxh = s_val - q2_val - aa1/2
    ok &= (ssh > 0) & (xxh > 0)

    return ok


# ── Method A: Full MC simulation ─────────────────────────────────────────
print("Running full MC simulation (sampling BW for each event)...")
rng = np.random.default_rng(42)
n_samples_per_evt = 10  # sample 10 masses per event to build statistics
mc_masses = []
n_report = N // 10

for j in range(N):
    if (j+1) % n_report == 0:
        print(f"  Event {j+1}/{N}...")

    am_max = ammax_per_evt[j]
    if am_max <= ammin:
        continue

    masses = sample_bw_python(am_max, n_samples_per_evt, rng)
    if len(masses) == 0:
        continue

    passed = check_all_fortran(masses, Q2[j], W2[j], t[j], sx[j], sxp[j], aly[j], s[j], x[j])
    mc_masses.extend(masses[passed].tolist())

mc_masses = np.array(mc_masses)
print(f"  MC simulation: {len(mc_masses)} accepted mass values from {N} events × {n_samples_per_evt} samples")

# ── Method B: Empirical BW fit ───────────────────────────────────────────
from scipy.optimize import curve_fit

def bw_fit_func(M, M0_fit, G0_fit, norm):
    """Running-width BW with free mass, width, and normalization."""
    M2 = M**2
    pm = np.sqrt(np.maximum(0, M2/4 - mpi**2))
    p0f = np.sqrt(max(0, M0_fit**2/4 - mpi**2))
    if p0f <= 0:
        return np.zeros_like(M)
    Grun = G0_fit * (pm/p0f)**3 * (M0_fit/M)
    return norm * 2*M * Grun / ((M2 - M0_fit**2)**2 + M0_fit**2 * Grun**2)


# ── Method C: BW × acceptance from event-by-event binary test ────────────
# For each M_test value, check what fraction of events would accept that M
def bw_dNdM(M):
    M2 = M**2
    pm = np.sqrt(np.maximum(0, M2/4 - mpi**2))
    Grun = G0 * (pm/p0)**3 * (M0/M)
    return 2*M * Grun / ((M2 - M0**2)**2 + M0**2 * Grun**2)

M_test = np.linspace(0.28, 1.10, 300)

# Probability that a random event can produce mass M:
# P(M) = BW(M) × <1(ammax_j > M) × 1(kin_checks_pass(M, event_j))> / <Z_j^{-1} × 1(ammax_j > M)>
# Actually the correct formula is:
# The effective dN/dM at mass M is proportional to:
#   BW(M) × sum_j [1(ammax_j > M) × 1(all_checks_pass)] / Z_j
# where Z_j = integral of BW from ammin to ammax_j

def bw_norm(ammax_val):
    if ammax_val <= ammin:
        return 0.0
    M_int = np.linspace(ammin, ammax_val, 1000)
    return np.trapz(bw_dNdM(M_int), M_int)

print("Computing per-event BW normalizations...")
unique_ammax = np.unique(np.round(ammax_per_evt, 5))
norm_lookup = {}
for am in unique_ammax:
    norm_lookup[am] = bw_norm(am)
Z_per_evt = np.array([norm_lookup[round(am, 5)] for am in ammax_per_evt])
Z_per_evt = np.where(Z_per_evt > 0, Z_per_evt, 1e-30)

eff_bw_Zweight = np.zeros_like(M_test)
eff_bw_noZ = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    can_produce = ammax_per_evt >= M_val
    acc = check_all_fortran(
        np.full(N, M_val), Q2, W2, t, sx, sxp, aly, s, x) & can_produce
    eff_bw_Zweight[i] = bw_dNdM(M_val) * (acc.astype(float) / Z_per_evt).sum()
    eff_bw_noZ[i] = bw_dNdM(M_val) * acc.sum() / N


# ── Plotting ──────────────────────────────────────────────────────────────
mrange = (0.28, 1.10)
nbins = 80
bw_width = (mrange[1]-mrange[0])/nbins
counts_data, edges = np.histogram(Mmeson, bins=nbins, range=mrange)
centers = 0.5*(edges[:-1] + edges[1:])
N_in = counts_data.sum()

# Also histogram the MC simulation
counts_mc, _ = np.histogram(mc_masses, bins=nbins, range=mrange)
N_mc_in = counts_mc.sum()

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Harut BW shape: MC simulation vs analytical models', fontsize=13)

def chi2_calc(data, model):
    mask = data > 0
    return np.sum((data[mask] - model[mask])**2 / data[mask])

# Panel 1: Data histogram + bare BW
ax = axes[0,0]
ax.bar(centers, counts_data, width=bw_width, color='lightskyblue', edgecolor='navy', alpha=0.7, label='LUND data')
bw_line = bw_dNdM(M_test)
integral_bw = np.trapz(bw_line, M_test)
sc_bw = N_in * bw_width / integral_bw
ax.plot(M_test, sc_bw*bw_line, '-', color='forestgreen', lw=2, label='Bare BW')
model_vals = np.interp(centers, M_test, sc_bw*bw_line)
c2 = chi2_calc(counts_data, model_vals)
ax.set_title(f'Bare BW  (χ²/ndf = {c2:.1f}/{nbins})')
ax.legend(fontsize=9)

# Panel 2: MC simulation overlay
ax = axes[0,1]
ax.bar(centers, counts_data, width=bw_width, color='lightskyblue', edgecolor='navy', alpha=0.7, label='LUND data')
if N_mc_in > 0:
    sc_mc = N_in / N_mc_in
    ax.bar(centers, sc_mc*counts_mc, width=bw_width*0.7, color='red', alpha=0.4, label=f'MC simulation (×{n_samples_per_evt}/evt)')
    c2_mc = chi2_calc(counts_data, sc_mc*counts_mc)
    ax.set_title(f'MC simulation  (χ²/ndf = {c2_mc:.1f}/{nbins})')
else:
    ax.set_title('MC simulation (no events)')
ax.legend(fontsize=9)

# Panel 3: BW × acceptance (no Z-weight)
ax = axes[1,0]
ax.bar(centers, counts_data, width=bw_width, color='lightskyblue', edgecolor='navy', alpha=0.7, label='LUND data')
integral_noZ = np.trapz(eff_bw_noZ, M_test)
if integral_noZ > 0:
    sc_noZ = N_in * bw_width / integral_noZ
    ax.plot(M_test, sc_noZ*eff_bw_noZ, '-', color='red', lw=2, label='BW × accept(M)')
    model_noZ = np.interp(centers, M_test, sc_noZ*eff_bw_noZ)
    c2_noZ = chi2_calc(counts_data, model_noZ)
    ax.set_title(f'BW × accept(M)  (χ²/ndf = {c2_noZ:.1f}/{nbins})')
ax.legend(fontsize=9)

# Panel 4: BW × accept / Z_j
ax = axes[1,1]
ax.bar(centers, counts_data, width=bw_width, color='lightskyblue', edgecolor='navy', alpha=0.7, label='LUND data')
integral_Z = np.trapz(eff_bw_Zweight, M_test)
if integral_Z > 0:
    sc_Z = N_in * bw_width / integral_Z
    ax.plot(M_test, sc_Z*eff_bw_Zweight, '-', color='purple', lw=2, label='BW × accept/Z_j')
    model_Z = np.interp(centers, M_test, sc_Z*eff_bw_Zweight)
    c2_Z = chi2_calc(counts_data, model_Z)
    ax.set_title(f'BW × accept/Z_j  (χ²/ndf = {c2_Z:.1f}/{nbins})')
ax.legend(fontsize=9)

for ax in axes.flat:
    ax.set_xlabel(r'$M_{\pi^+\pi^-}$ (GeV)')
    ax.set_ylabel('Events')
    ax.grid(True, alpha=0.4, ls='--')

fig.tight_layout()
fig.savefig('/tmp/harut_bw_mc.png', dpi=150)
print("Saved /tmp/harut_bw_mc.png")

# ── Also try a BW fit with free parameters ──────────────────────────────
fig2, ax2 = plt.subplots(figsize=(9, 6))
mfit_range = (0.50, 1.05)
nbins2 = 80
bw_width2 = (mfit_range[1]-mfit_range[0])/nbins2
counts2, edges2 = np.histogram(Mmeson, bins=nbins2, range=mfit_range)
centers2 = 0.5*(edges2[:-1] + edges2[1:])
N_in2 = counts2.sum()

ax2.bar(centers2, counts2, width=bw_width2, color='lightskyblue', edgecolor='navy', alpha=0.7, label='LUND data')

# Fit running-width BW with free M0, G0, norm
try:
    popt, pcov = curve_fit(bw_fit_func, centers2, counts2,
                           p0=[M0, G0, N_in2*bw_width2],
                           bounds=([0.5, 0.05, 0], [1.0, 0.5, N_in2*100]),
                           sigma=np.sqrt(np.where(counts2 > 0, counts2, 1)))
    M0_fit, G0_fit, norm_fit = popt
    print(f"\nBW fit: M0 = {M0_fit:.5f} GeV, G0 = {G0_fit:.5f} GeV")

    M_fit_line = np.linspace(mfit_range[0], mfit_range[1], 500)
    ax2.plot(M_fit_line, bw_fit_func(M_fit_line, *popt), '-', color='red', lw=2,
             label=f'BW fit: M₀={M0_fit:.4f}, Γ₀={G0_fit:.4f}')

    model_fit = bw_fit_func(centers2, *popt)
    mask = counts2 > 0
    c2_fit = np.sum((counts2[mask] - model_fit[mask])**2 / counts2[mask])

    # Also overlay nominal BW
    bw_nom = bw_dNdM(M_fit_line)
    integral_nom = np.trapz(bw_nom, M_fit_line)
    sc_nom = N_in2 * bw_width2 / integral_nom
    ax2.plot(M_fit_line, sc_nom*bw_nom, '--', color='forestgreen', lw=1.5, alpha=0.7,
             label=f'Nominal BW: M₀={M0:.5f}, Γ₀={G0:.4f}')

    ax2.set_title(f'BW fit  (χ²/ndf = {c2_fit:.1f}/{nbins2}, M₀={M0_fit:.4f}, Γ₀={G0_fit:.4f})')
    print(f"  Fit chi2/ndf = {c2_fit:.1f}/{nbins2}")
except Exception as e:
    print(f"Fit failed: {e}")
    ax2.set_title('BW fit failed')

ax2.set_xlabel(r'$M_{\pi^+\pi^-}$ (GeV)', fontsize=13)
ax2.set_ylabel('Events')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.4, ls='--')
fig2.tight_layout()
fig2.savefig('/tmp/harut_bw_fit.png', dpi=150)
print("Saved /tmp/harut_bw_fit.png")

import os
os.system('open /tmp/harut_bw_mc.png /tmp/harut_bw_fit.png')
