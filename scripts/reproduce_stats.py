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
from scipy.optimize import brentq

# Pre-specified equivalence margin (two one-sided tests). Choose on biological
# grounds -- the largest difference that would be considered negligible -- NOT to
# obtain significance. Expressed as a percentage of the 2.2 kb reference mean.
EQUIV_MARGIN_PCT = 20.0

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

# --- Equivalence (TOST / 90% CI of the difference) and power ---
# A non-significant difference is not evidence of equivalence. Two one-sided tests
# (TOST) at alpha=0.05 are equivalent to checking that the 90% CI of the difference
# lies entirely within the pre-specified margin.
sp = np.sqrt(((len(a) - 1) * a.std(ddof=1) ** 2 + (len(b) - 1) * b.std(ddof=1) ** 2)
             / (len(a) + len(b) - 2))
se = sp * np.sqrt(1 / len(a) + 1 / len(b))
dfree = len(a) + len(b) - 2
diff = b.mean() - a.mean()
ref = a.mean()                                   # express relative to the 2.2 kb mean
tcrit90 = stats.t.ppf(0.95, dfree)               # 90% CI <-> TOST at alpha = 0.05
ci_lo, ci_hi = diff - tcrit90 * se, diff + tcrit90 * se
margin = EQUIV_MARGIN_PCT / 100 * ref
p_tost = max(stats.t.sf((diff + margin) / se, dfree), stats.t.cdf((diff - margin) / se, dfree))


def _power(dd):
    ncp = abs(dd) * np.sqrt(len(a) * len(b) / (len(a) + len(b)))
    tc = stats.t.ppf(0.975, dfree)
    return float(1 - stats.nct.cdf(tc, dfree, ncp) + stats.nct.cdf(-tc, dfree, ncp))


mde = brentq(lambda x: _power(x) - 0.80, 0.05, 1.4)
print("\nequivalence / power:")
print(f"  difference (2.6-2.2) = {diff:+.1f} A.U. ({diff / ref * 100:+.1f}%)")
print(f"  90% CI of difference = [{ci_lo:+.1f}, {ci_hi:+.1f}] A.U. "
      f"= [{ci_lo / ref * 100:+.1f}%, {ci_hi / ref * 100:+.1f}%]")
print(f"  TOST equivalence within +/-{EQUIV_MARGIN_PCT:.0f}%: p = {p_tost:.3f} "
      f"-> {'EQUIVALENT' if p_tost < 0.05 else 'not shown'}")
print(f"  achieved power (observed d={d:.2f}) = {_power(d) * 100:.0f}%; "
      f"min. detectable difference at 80% power = {mde * sp:.0f} A.U. ({mde * sp / ref * 100:.0f}%)")

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
