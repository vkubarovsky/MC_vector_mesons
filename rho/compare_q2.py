"""
compare_q2.py — overlay Q2, y, W distributions for DIFFRAD vs Harut sigma_T
with narrow xB and t cuts, to isolate the Q2 dependence of sigma_T.

Cuts: Q2=[1,5.5], xB=[0.24,0.26], -t=[0.4,0.6], y=[0.10,0.80], W>=2.0
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lund_reader import read_lund_events

DIFFRAD_LUND = os.path.expanduser('~/Downloads/DIFFRAD_lund/rho/born_events.lund')
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

print("Reading DIFFRAD ...", flush=True)
ev_d = read_lund_events(DIFFRAD_LUND)
D = apply_cuts(ev_d, cuts)
print(f"  {len(ev_d)} total → {len(D['Q2'])} pass cuts")

print("Reading Harut sigma_T ...", flush=True)
ev_h = read_lund_events(HARUT_LUND)
H = apply_cuts(ev_h, cuts)
print(f"  {len(ev_h)} total → {len(H['Q2'])} pass cuts")

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
    nd, nh = len(D[key]), len(H[key])
    ax.hist(D[key], bins=bins, color='tab:blue',
            label=f"DIFFRAD (E=10.6, N={nd})", **kw)
    ax.hist(H[key], bins=bins, color='tab:green',
            label=f"Harut σ_T (E=10.2, N={nh})", **kw)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Normalized', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)
    # annotate means
    ax.axvline(D[key].mean(), color='tab:blue',  ls='--', lw=1,
               label=f'⟨D⟩={D[key].mean():.3f}')
    ax.axvline(H[key].mean(), color='tab:green', ls='--', lw=1,
               label=f'⟨H⟩={H[key].mean():.3f}')
    ax.legend(fontsize=8)

fig.suptitle(f'DIFFRAD vs Harut σ_T  |  Cuts: {cut_str}', fontsize=10)
fig.tight_layout()

out = os.path.join(os.path.dirname(__file__), 'compare_q2.png')
fig.savefig(out, dpi=150)
print(f"\nSaved: {out}")

for key in ('Q2', 'y', 'W'):
    print(f"  {key}: DIFFRAD ⟨{key}⟩={D[key].mean():.4f}  Harut ⟨{key}⟩={H[key].mean():.4f}")
