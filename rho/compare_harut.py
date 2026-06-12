"""
compare_harut.py — overlay kinematic distributions from three generators:
  1. DIFFRAD (born)
  2. gagrho (fixed, dsigmaL=0) — full and with DIFFRAD cuts
  3. Harut's generator

Saves: rho/harut_comparison.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lund_reader import read_lund_events

DIFFRAD_LUND = os.path.expanduser('~/Downloads/DIFFRAD_lund/rho/born_events.lund')
GAGRHO_LUND  = os.path.expanduser('~/Downloads/rho_gagik_fixed.lund')
HARUT_LUND   = '/Users/vpk/Downloads/Harut_rho_MC_generator/Harut_sigma_U.lund'

print("Reading DIFFRAD ...", flush=True)
ev_d = read_lund_events(DIFFRAD_LUND)
print("Reading gagrho (fixed) ...", flush=True)
ev_g = read_lund_events(GAGRHO_LUND)
print("Reading Harut ...", flush=True)
ev_h = read_lund_events(HARUT_LUND)

def arr(evs, key):
    return np.array([e[key] for e in evs])

def extract(evs):
    return {k: arr(evs, k) for k in ('Q2', 'xB', 'W', 'y', 't')}

D = extract(ev_d)
G = extract(ev_g)
H = extract(ev_h)

# Apply DIFFRAD kinematic cuts to gagrho
mask = ((G['Q2'] <= D['Q2'].max()) & (G['W'] >= D['W'].min()) &
        (G['xB'] >= D['xB'].min()) & (G['xB'] <= D['xB'].max()))
G_cut = {k: G[k][mask] for k in G}

datasets = [
    (D,     'tab:blue',   f"DIFFRAD  (E=10.6, N={len(ev_d)//1000}k, ⟨y⟩={D['y'].mean():.3f})"),
    (G,     'tab:orange', f"gagrho   (E=10.6, N={len(ev_g)//1000}k, ⟨y⟩={G['y'].mean():.3f})"),
    (G_cut, 'tab:red',    f"gagrho+cuts (⟨y⟩={G_cut['y'].mean():.3f})"),
    (H,     'tab:green',  f"Harut    (E=10.2, N={len(ev_h)//1000}k, ⟨y⟩={H['y'].mean():.3f})"),
]

vars_info = [
    ('Q2', r'$Q^2$ (GeV²)',  (0, 10),   50),
    ('xB', r'$x_B$',         (0, 0.75), 50),
    ('W',  r'$W$ (GeV)',     (1.5, 5),  50),
    ('y',  r'$y$',           (0, 1),    50),
    ('t',  r'$-t$ (GeV²)',   (0, 6),    50),
]

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
axes = axes.flatten()

for ax, (key, xlabel, xlim, nbins) in zip(axes, vars_info):
    bins = np.linspace(xlim[0], xlim[1], nbins + 1)
    for d, col, label in datasets:
        vals = -d['t'] if key == 't' else d[key]
        ls = '--' if 'cuts' in label else '-'
        ax.hist(vals, bins=bins, histtype='step', density=True,
                color=col, linewidth=1.5, linestyle=ls, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Normalized events')
    ax.set_xlim(xlim)
    ax.legend(fontsize=7)

# 6th panel: y gap decomposition as bar chart
ax = axes[5]
labels = ['gagrho\nfull', 'gagrho\n+W,Q² cuts', 'DIFFRAD']
means  = [G['y'].mean(), G_cut['y'].mean(), D['y'].mean()]
colors = ['tab:orange', 'tab:red', 'tab:blue']
bars = ax.bar(labels, means, color=colors, width=0.5)
ax.set_ylabel(r'$\langle y \rangle$')
ax.set_ylim(0.35, 0.52)
ax.set_title(r'$\langle y \rangle$ gap decomposition')
for bar, val in zip(bars, means):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.003, f'{val:.3f}',
            ha='center', va='bottom', fontsize=10, fontweight='bold')
# annotate arrows
dy_cuts = G_cut['y'].mean() - G['y'].mean()
dy_model = D['y'].mean() - G_cut['y'].mean()
ax.annotate(f'W_min cut\n+{dy_cuts:.3f}', xy=(1, G['y'].mean()+dy_cuts/2),
            xytext=(1.5, G['y'].mean()+dy_cuts/2), fontsize=8, color='tab:red',
            arrowprops=dict(arrowstyle='->', color='tab:red'))
ax.annotate(f'σ model\n+{dy_model:.3f}', xy=(1.5, G_cut['y'].mean()+dy_model/2),
            xytext=(2.1, G_cut['y'].mean()+dy_model/2+0.01), fontsize=8, color='tab:blue',
            arrowprops=dict(arrowstyle='->', color='tab:blue'))

fig.suptitle('ρ kinematic comparison — DIFFRAD vs gagrho (fixed) vs Harut\n'
             f'DIFFRAD cuts: Q²<{D["Q2"].max():.1f}, W>{D["W"].min():.2f} GeV  '
             f'(gap: {D["y"].mean()-G["y"].mean():.3f} total = {dy_cuts:.3f} kin + {dy_model:.3f} model)',
             fontsize=11)
fig.tight_layout()

outfile = os.path.join(os.path.dirname(__file__), 'harut_comparison.png')
fig.savefig(outfile, dpi=150)
print(f"\nSaved: {outfile}")

print(f"\n{'='*60}")
print(f"<y> gap decomposition:")
print(f"  Total gap (DIFFRAD - gagrho):     {D['y'].mean()-G['y'].mean():.4f}")
print(f"  From kin cuts (W_min, Q2_max):    {dy_cuts:.4f}  ({dy_cuts/(D['y'].mean()-G['y'].mean())*100:.0f}%)")
print(f"  From cross-section model:         {dy_model:.4f}  ({dy_model/(D['y'].mean()-G['y'].mean())*100:.0f}%)")
print(f"\nAt same Q², DIFFRAD has lower <xB> → higher <y> (from difflt σ_T model)")
