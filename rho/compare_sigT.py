"""
compare_sigT.py — three-way Q² comparison to isolate sigma_L contamination.

Samples compared (all E=10.2 GeV, same narrow cuts):
  1. DIFFRAD born (sigma_T + eps*sigma_L) — ~/Downloads/DIFFRAD_lund/rho/born_events.lund
  2. DIFFRAD sigT=0 (sigma_T only, sigmal=0 in bornin) — sigT_born_events.lund  [local]
  3. Harut sigma_T — ~/Downloads/Harut_rho_MC_generator/Harut_sigma_T.lund

Cuts: Q2=[1,5.5], xB=[0.24,0.26], -t=[0.4,0.6], y=[0.10,0.80], W>=2.0
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    return dict(Q2=Q2[mask], xB=xB[mask], y=y[mask], W=W[mask], t=t[mask])

print("Reading DIFFRAD (sigT+sigL) ...", flush=True)
ev_d = read_lund_events(DIFFRAD_LUND)
D = apply_cuts(ev_d, cuts)
print(f"  {len(ev_d)} total → {len(D['Q2'])} pass cuts")

print("Reading DIFFRAD (sigT only) ...", flush=True)
ev_s = read_lund_events(SIGT_LUND)
S = apply_cuts(ev_s, cuts)
print(f"  {len(ev_s)} total → {len(S['Q2'])} pass cuts")

print("Reading Harut sigma_T ...", flush=True)
ev_h = read_lund_events(HARUT_LUND)
H = apply_cuts(ev_h, cuts)
print(f"  {len(ev_h)} total → {len(H['Q2'])} pass cuts")

# Theoretical <Q2> for pure sigma_T (flat in Q2 at fixed xB → weight ~ 1/(1-eps))
# Computed numerically for xB=0.25, E=10.2 GeV
Mp = 0.938272
E  = 10.2
s  = 2*Mp*E
xB_mid = 0.25
q2_vals = np.linspace(1.0, 5.5, 5000)
y_vals  = q2_vals / (s * xB_mid)
valid   = (y_vals > 0.1) & (y_vals < 0.8)
q2_v    = q2_vals[valid]; y_v = y_vals[valid]
eps_v   = (1-y_v) / (1 - y_v + 0.5*y_v**2)
# Weight for sigma_T only: dN/dQ2 ~ 1/(1-eps)
wt      = 1.0 / (1.0 - eps_v)
q2_theory = np.average(q2_v, weights=wt)
print(f"\nTheoretical <Q2> (pure sigma_T, xB=0.25, E=10.2): {q2_theory:.3f} GeV2")

datasets = [
    (D, 'tab:blue',   f"DIFFRAD σT+εσL (N={len(D['Q2'])})"),
    (S, 'tab:orange', f"DIFFRAD σT only (N={len(S['Q2'])})"),
    (H, 'tab:green',  f"Harut σT       (N={len(H['Q2'])})"),
]

vars_info = [
    ('Q2', r'$Q^2$ (GeV²)', (1.0, 5.5), 20),
    ('y',  r'$y$',          (0.1, 0.8), 20),
    ('W',  r'$W$ (GeV)',    (2.0, 4.0), 20),
    ('xB', r'$x_B$',        (0.23, 0.27), 10),
]

fig, axes = plt.subplots(1, 4, figsize=(18, 4))
cut_str = (f"xB∈[{cuts['xbmin']},{cuts['xbmax']}]  "
           f"−t∈[{cuts['tmin']},{cuts['tmax']}]  "
           f"Q²∈[{cuts['q2min']},{cuts['q2max']}]  "
           f"y∈[{cuts['ymin']},{cuts['ymax']}]  W≥{cuts['wmin']}")

for ax, (key, xlabel, xlim, nbins) in zip(axes, vars_info):
    bins = np.linspace(xlim[0], xlim[1], nbins + 1)
    kw   = dict(histtype='step', density=True, linewidth=2)
    for d, col, label in datasets:
        ax.hist(d[key], bins=bins, color=col, label=label, **kw)
    # vertical mean lines
    for d, col, _ in datasets:
        ax.axvline(d[key].mean(), color=col, ls='--', lw=1.2)
    if key == 'Q2':
        ax.axvline(q2_theory, color='black', ls=':', lw=1.5,
                   label=f'Theory pure σT: {q2_theory:.3f}')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Normalized', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

# Annotate means as text on Q2 panel
ax0 = axes[0]
for i, (d, col, label) in enumerate(datasets):
    ax0.text(0.97, 0.97 - i*0.12,
             f'⟨Q²⟩={d["Q2"].mean():.3f}', color=col, fontsize=9,
             ha='right', va='top', transform=ax0.transAxes)
ax0.text(0.97, 0.97 - 3*0.12,
         f'Theory={q2_theory:.3f}', color='black', fontsize=9,
         ha='right', va='top', transform=ax0.transAxes)

fig.suptitle(f'σT vs σT+εσL  |  E=10.2 GeV  |  Cuts: {cut_str}', fontsize=10)
fig.tight_layout()

out = os.path.join(HERE, 'compare_sigT.png')
fig.savefig(out, dpi=150)
print(f"\nSaved: {out}")

print(f"\n{'='*60}")
print(f"{'Sample':<30}  {'<Q2>':<8}  {'<y>':<8}  {'<W>':<8}  N")
print(f"{'-'*60}")
labels_short = ['DIFFRAD σT+εσL', 'DIFFRAD σT only', 'Harut σT']
for (d, col, _), lbl in zip(datasets, labels_short):
    print(f"  {lbl:<28}  {d['Q2'].mean():.4f}   {d['y'].mean():.4f}   {d['W'].mean():.4f}  {len(d['Q2'])}")
print(f"  {'Theory (pure σT)':<28}  {q2_theory:.4f}")
