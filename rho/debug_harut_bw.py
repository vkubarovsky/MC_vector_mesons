#!/usr/bin/env python3
"""
debug_harut_bw.py -- Compute the effective BW shape for Harut events by
replicating the exact accept/reject chain from diffrad_harut.f90:
  1. BW sampling with per-event ammax
  2. Threshold: W > Mp + M
  3. Discriminant checks (disc1, disc2 >= 0)
  4. t in kinematic range [tdmink, tdmaxk]
  5. vmax > 0 check
  6. ssh > 0 and xxh > 0 check
"""

import numpy as np
import matplotlib.pyplot as plt
from lund_reader import read_lund_events

# BW parameters
M0   = 0.77526
G0   = 0.1502
mpi  = 0.13957
Mp   = 0.938272
Mp2  = Mp**2
ammin = 2*mpi

p0 = np.sqrt(max(0, M0**2/4 - mpi**2))

# Read events
fname = '/Users/vpk/Downloads/DIFFRAD_lund/harut/rho/born_events.lund'
print(f"Reading {fname} ...")
events = read_lund_events(fname, max_ev=50000)
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

# Derived Fortran variables (from conkin)
s   = 2*Mp*ebeam + Mp2
ys  = Q2 / (s * xB)
x   = s * (1 - ys)
sx  = s - x
sxp = s + x
aly = sx**2 + 4*Mp2*Q2

# tdif in Fortran is the Mandelstam t (negative); in our LUND t = 2*Mp*(Mp - Ep)
# which is t = (p'-p)^2 < 0. The Fortran tdif is sampled and should equal t.
tdif = t

# Per-event BW upper limit
ammax_per_evt = np.minimum(W - Mp, M0 + 5*G0)

print(f"  ammax range: {ammax_per_evt.min():.4f} – {ammax_per_evt.max():.4f} GeV")
print(f"  M0+5*G0 = {M0+5*G0:.4f} GeV")
print(f"  Fraction with ammax < M0+5*G0: {(ammax_per_evt < M0+5*G0 - 0.001).sum()/N*100:.1f}%")

# Running-width BW: dN/dM
def bw_dNdM(M):
    M2 = M**2
    pm = np.sqrt(np.maximum(0, M2/4 - mpi**2))
    Grun = G0 * (pm/p0)**3 * (M0/M)
    return 2*M * Grun / ((M2 - M0**2)**2 + M0**2 * Grun**2)

# BW normalization integral for a given ammax
def bw_norm(ammax_val):
    if ammax_val <= ammin:
        return 0.0
    M_int = np.linspace(ammin, ammax_val, 500)
    return np.trapz(bw_dNdM(M_int), M_int)

# Full kinematic acceptance including all M-dependent checks from diffrad_harut
def full_accept(M_val, q2, w2, t_arr, sx_arr, sxp_arr, aly_arr, s_arr, x_arr):
    """Check all M-dependent filters from diffrad_harut.f90 main loop."""
    amv2 = M_val**2

    # Check 1: W > Mp + M  (line 328)
    ok = (np.sqrt(w2) >= Mp + M_val)

    # Check 2: disc1, disc2 >= 0  (lines 334-337)
    tt1 = w2 - q2 - Mp2
    tt2 = w2 - Mp2 + amv2
    disc1 = tt1**2 + 4*q2*w2
    disc2 = tt2**2 - 4*amv2*w2
    ok &= (disc1 >= 0) & (disc2 >= 0)

    # Check 3: t in [tdmink, tdmaxk]  (lines 339-342)
    safe = ok.copy()
    sd1 = np.sqrt(np.maximum(0, disc1))
    sd2 = np.sqrt(np.maximum(0, disc2))
    tdmink = np.where(safe, -q2 + amv2 - 0.5/w2*(tt1*tt2 + sd1*sd2), 0)
    tdmaxk = np.where(safe, -q2 + amv2 - 0.5/w2*(tt1*tt2 - sd1*sd2), 0)
    ok &= (t_arr >= tdmink) & (t_arr <= tdmaxk)

    # Check 4: vmax > 0  (lines 346-363)
    sxt = sx_arr + t_arr
    tq  = q2 + t_arr - amv2
    vmax_kin = np.where(ok,
        tt2 + 0.5/q2*(-tt1*tq
            - np.sqrt(np.maximum(0, tt1**2 + 4*q2*w2))
            * np.sqrt(np.maximum(0, tq**2 + 4*amv2*q2)))
        - 1e-8,
        1.0)  # dummy positive for non-ok events
    ok &= (vmax_kin > 0)

    # Check 5: ssh > 0 and xxh > 0  (lines 373-378)
    aa1 = np.where(ok,
        (q2*sxp_arr*sxt - (s_arr*sx_arr + 2*Mp2*q2)*tq) / (2*aly_arr), 0)
    aa2 = np.where(ok,
        (q2*sxp_arr*sxt - (x_arr*sx_arr - 2*Mp2*q2)*tq) / (2*aly_arr), 0)
    vv1 = aa1 / 2
    vv2 = aa2 / 2
    ssh = x_arr + q2 - vv2
    xxh = s_arr - q2 - vv1
    ok &= (ssh > 0) & (xxh > 0)

    return ok


M_test = np.linspace(0.60, 0.95, 200)

# ── Method 1: Bare BW ────────────────────────────────────────────────────
bw_bare = bw_dNdM(M_test)

# ── Method 2: BW × basic kinematic acceptance (checks 1-3) ──────────────
def basic_kin_accept(M_val, q2, w2, t_arr):
    amv2 = M_val**2
    ok = (np.sqrt(w2) >= Mp + M_val)
    tt1 = w2 - q2 - Mp2
    tt2 = w2 - Mp2 + amv2
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

eff_bw_basic = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    acc = basic_kin_accept(M_val, Q2, W2, tdif)
    eff_bw_basic[i] = bw_dNdM(M_val) * acc.sum() / N

# ── Method 3: BW × ALL acceptance checks (1-5) ──────────────────────────
eff_bw_allcuts = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    acc = full_accept(M_val, Q2, W2, tdif, sx, sxp, aly, s, x)
    eff_bw_allcuts[i] = bw_dNdM(M_val) * acc.sum() / N

# ── Method 4: BW × ALL checks × per-event 1/Z_j normalization ───────────
print("Computing per-event BW normalizations...")
unique_ammax = np.unique(np.round(ammax_per_evt, 5))
norm_lookup = {}
for am in unique_ammax:
    norm_lookup[am] = bw_norm(am)

Z_per_evt = np.array([norm_lookup[round(am, 5)] for am in ammax_per_evt])
Z_per_evt = np.where(Z_per_evt > 0, Z_per_evt, 1e-30)
Z_full = bw_norm(M0 + 5*G0)

eff_bw_full = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    can_produce = ammax_per_evt > M_val
    acc = full_accept(M_val, Q2, W2, tdif, sx, sxp, aly, s, x) & can_produce
    weights = acc.astype(float) / Z_per_evt
    eff_bw_full[i] = bw_dNdM(M_val) * weights.sum()

# ── Plot comparison ──────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Harut BW shape analysis: all M-dependent checks from diffrad_harut.f90', fontsize=13)

mrange = (0.60, 0.95)
nbins = 60
bw_width = (mrange[1]-mrange[0])/nbins
counts, edges = np.histogram(Mmeson, bins=nbins, range=mrange)
centers = 0.5*(edges[:-1] + edges[1:])
N_in = counts.sum()

def plot_overlay(ax, M_arr, y_arr, color, label, title_prefix):
    ax.bar(centers, counts, width=bw_width, color='lightskyblue', edgecolor='navy', alpha=0.7)
    integral = np.trapz(y_arr, M_arr)
    if integral > 0:
        sc = N_in * bw_width / integral
        ax.plot(M_arr, sc*y_arr, '-', color=color, lw=2, label=label)
        chi2 = np.sum((counts - np.interp(centers, M_arr, sc*y_arr))**2 /
                      np.where(counts > 0, counts, 1))
        ax.set_title(f'{title_prefix}  (χ²/ndf = {chi2:.1f}/{nbins})')
        print(f"  {title_prefix}: chi2/ndf = {chi2:.1f}/{nbins}")
    else:
        ax.set_title(f'{title_prefix}  (failed)')
    ax.legend()

plot_overlay(axes[0,0], M_test, bw_bare, 'forestgreen', 'Bare BW', 'Bare BW')
plot_overlay(axes[0,1], M_test, eff_bw_basic, 'red', 'BW × basic kin (1-3)', 'BW × basic kin')
plot_overlay(axes[1,0], M_test, eff_bw_allcuts, 'darkred', 'BW × all cuts (1-5)', 'BW × all cuts')
plot_overlay(axes[1,1], M_test, eff_bw_full, 'purple', 'BW × all cuts × 1/Z_j', 'BW × all cuts × 1/Z_j')

for ax in axes.flat:
    ax.set_xlabel(r'$M_{\pi^+\pi^-}$ (GeV)')
    ax.set_ylabel('Events')
    ax.grid(True, alpha=0.4, ls='--')

fig.tight_layout()
fig.savefig('/tmp/harut_bw_methods.png', dpi=150)
print("Saved /tmp/harut_bw_methods.png")

# ── Acceptance fraction comparison ──────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))
fig2.suptitle('Acceptance fractions vs M', fontsize=12)

acc_basic = np.zeros_like(M_test)
acc_full  = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    acc_basic[i] = basic_kin_accept(M_val, Q2, W2, tdif).sum() / N
    acc_full[i]  = full_accept(M_val, Q2, W2, tdif, sx, sxp, aly, s, x).sum() / N

axes2[0].plot(M_test, acc_basic, '-b', lw=2, label='Basic kin (1-3)')
axes2[0].plot(M_test, acc_full, '-r', lw=2, label='All cuts (1-5)')
axes2[0].set_xlabel('M (GeV)'); axes2[0].set_ylabel('Acceptance fraction')
axes2[0].set_title('Acceptance vs M'); axes2[0].legend()

# Ratio: how much do vmax+ssh/xxh cuts add?
ratio = np.where(acc_basic > 0, acc_full / acc_basic, 1.0)
axes2[1].plot(M_test, ratio, '-k', lw=2)
axes2[1].set_xlabel('M (GeV)'); axes2[1].set_ylabel('Full / Basic')
axes2[1].set_title('Extra rejection from vmax + ssh/xxh')
axes2[1].axhline(1.0, color='r', ls='--')

# Fraction excluded by BW upper limit
ammax_frac = np.zeros_like(M_test)
for i, M_val in enumerate(M_test):
    ammax_frac[i] = (ammax_per_evt < M_val).sum() / N
axes2[2].plot(M_test, ammax_frac, '-r', lw=2)
axes2[2].set_xlabel('M (GeV)'); axes2[2].set_ylabel('Fraction with ammax < M')
axes2[2].set_title('Events excluded by BW upper limit')

for ax in axes2:
    ax.grid(True, alpha=0.4, ls='--')

fig2.tight_layout()
fig2.savefig('/tmp/harut_bw_corrections.png', dpi=150)
print("Saved /tmp/harut_bw_corrections.png")

import os
os.system('open /tmp/harut_bw_methods.png /tmp/harut_bw_corrections.png')
