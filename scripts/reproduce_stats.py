"""Reproduce the quantification (statistics + violin plot) from the committed
per-fish measurement tables -- no raw image data required.

    python scripts/reproduce_stats.py

Reads measurements/measurements_14a-ef1a.csv (2.2 kb) and
measurements/measurements_14a-ef1a-extra.csv (2.6 kb), each one row per fish, and
regenerates the primary result: mean-over-stack intensity compared between the two
EF1A promoter lengths. Writes figures/quantification.png and .pdf.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams['pdf.fonttype'] = 42          # Illustrator-editable text
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parent.parent
MEAS = ROOT / 'measurements'
FIGS = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

PRIMARY = 'stack_mean'                      # mean over all voxels of the bg-subtracted stack
GROUPS = [('2.2 kb', 'measurements_14a-ef1a.csv'),
          ('2.6 kb', 'measurements_14a-ef1a-extra.csv')]

data, summary = {}, []
for label, fn in GROUPS:
    df = pd.read_csv(MEAS / fn)
    v = df[PRIMARY].to_numpy(dtype=float)
    data[label] = v
    summary.append(dict(promoter=label, n=len(v), mean=v.mean(), sd=v.std(ddof=1),
                        mip_mean=df['mip_mean'].mean()))

s = pd.DataFrame(summary).set_index('promoter').round(2)
print(s.to_string())

a, b = data['2.2 kb'], data['2.6 kb']
u, p = stats.mannwhitneyu(a, b, alternative='two-sided')
t, pt = stats.ttest_ind(a, b)
d = (b.mean() - a.mean()) / np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
print(f"\nmean-over-stack (primary metric), fish as the unit:")
print(f"  n = {len(a)} (2.2 kb) vs {len(b)} (2.6 kb)")
print(f"  fold (2.6/2.2) = {b.mean() / a.mean():.2f}x")
print(f"  Mann-Whitney U = {u:.1f}, p = {p:.3f}")
print(f"  Welch t-test   t = {t:.2f}, p = {pt:.3f}")
print(f"  Cohen's d      = {d:.2f}")
sig = 'n.s.' if p >= 0.05 else (f'p = {p:.3f}' if p >= 1e-3 else f'p = {p:.1e}')

# --- violin ---
order = list(data)
vals = [data[k] for k in order]
fig, ax = plt.subplots(figsize=(4.5, 5))
parts = ax.violinplot(vals, showextrema=False)
for pc in parts['bodies']:
    pc.set_facecolor('white'); pc.set_edgecolor('black'); pc.set_alpha(1.0); pc.set_linewidth(1.0)
rng = np.random.default_rng(0)
for i, g in enumerate(vals):
    ax.scatter(rng.normal(i + 1, 0.05, len(g)), g, s=22, color='gray', alpha=0.8, zorder=2)
    q1, med, q3 = np.percentile(g, [25, 50, 75])
    iqr = q3 - q1
    lo, hi = max(g.min(), q1 - 1.5 * iqr), min(g.max(), q3 + 1.5 * iqr)
    ax.plot([i + 1, i + 1], [lo, hi], color='black', lw=1.2, zorder=3)
    ax.plot([i + 1, i + 1], [q1, q3], color='black', lw=6, zorder=3)
    ax.scatter([i + 1], [med], s=90, facecolor='white', edgecolor='black', linewidth=1.5, zorder=4)
ytop = max(g.max() for g in vals)
ax.plot([1, 1, 2, 2], [ytop * 1.03, ytop * 1.06, ytop * 1.06, ytop * 1.03], color='black', lw=1.0)
ax.text(1.5, ytop * 1.065, sig, ha='center', va='bottom', fontsize=12)
ax.set_xticks([1, 2]); ax.set_xticklabels([f'{k}\n(n={len(data[k])})' for k in order], fontsize=12)
ax.set_ylabel('mean-over-stack intensity (A.U.)', fontsize=12)
ax.set_xlabel('EF1A promoter length', fontsize=12)
ax.set_ylim(0, ytop * 1.15)
ax.spines[['top', 'right']].set_visible(False)
ax.set_title('p14a.eef1a1l1:mSG:sec61b', fontsize=12)
plt.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(FIGS / f'quantification.{ext}', dpi=300, bbox_inches='tight')
print(f"\nWrote {FIGS / 'quantification.png'} and .pdf")
