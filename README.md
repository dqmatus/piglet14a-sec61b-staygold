# EF1A promoter-length comparison at pIGLET14a (StayGold–Sec61β)

Single-channel confocal quantification comparing two lengths of the zebrafish
**EF1A (`eef1a1l1`) promoter** — **2.2 kb** vs **2.6 kb** — driving an
ER-membrane reporter, **monomeric StayGold fused to Sec61β**, from the
**pIGLET14a** genomic landing site. Construct: `p14a.eef1a1l1:mSG:sec61b`.

**Result.** On uniformly acquired datasets with the *fish* as the unit of
replication (2.2 kb, n = 17; 2.6 kb, n = 30), the two promoter lengths drive
**statistically indistinguishable** reporter expression:

| promoter | n (fish) | mean-over-stack (A.U.) |
|----------|:--------:|:----------------------:|
| 2.2 kb   |    17    |       482 ± 49         |
| 2.6 kb   |    30    |       464 ± 37         |

Mann–Whitney **p = 0.26** (primary metric); Welch *t*-test p = 0.15; Cohen's
*d* = −0.43; fold (2.6/2.2) = 0.96×. The legacy max-intensity-projection metric
agrees (748 vs 715 A.U., p = 0.25). This is a well-powered null: the additional
0.4 kb of EF1A promoter does not measurably change expression at pIGLET14a.

## Constructs

| folder / group | promoter | fish |
|----------------|----------|------|
| `14a-ef1a`       | pIGLET14a EF1A **2.2 kb** | 17 (1 FOV each) |
| `14a-ef1a-extra` | pIGLET14a EF1A **2.6 kb** | 30 (1 FOV each) |

Every fish is one image (one field of view), imaged in the same anatomical
region (zebrafish flank, posterior to the yolk extension) at the same stage
(72 hpf) with identical acquisition settings. The fish — not the field of view —
is the unit of replication.

## Quantification

1. **Background subtraction** — a dark/light calibration model from laser-off
   (camera offset) and laser-on/no-sample (total background) reference
   acquisitions, subtracted per pixel and clipped at zero.
2. **Primary metric — `stack_mean` (mean-over-stack):** the mean over all voxels
   of the background-subtracted z-stack. It is linear (∝ total signal ÷ z-plane
   count), projection-independent, and insensitive to the number of z-planes.
3. **Robustness cross-check — `mip_mean`:** the mean of the max-intensity
   projection (the conventional metric). Reported alongside to confirm the result
   is not a projection artifact.
4. **Comparison** — per-fish values, Mann–Whitney U (primary), Welch *t*-test,
   fold change, Cohen's *d*, and a figure.

Single channel, **no internal control** — absolute intensities rely on matched
acquisition across sessions (confirmed identical here).

## Repository layout

```
notebooks/    Full analysis notebook (raw ND2 -> per-fish measurements -> figures)
measurements/ Per-fish and per-FOV measurement tables (committed)
figures/      Publication figure (PNG + Illustrator-editable PDF)
scripts/      reproduce_stats.py — regenerate stats + plot from the CSVs
```

## Reproducing the analysis

**From the committed measurement tables (no raw data needed)** — regenerates the
statistics and the quantification plot:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python scripts/reproduce_stats.py
```

**Full pipeline from raw images** — `notebooks/01_ef1a_sec61b_quantification.ipynb`
runs background subtraction, per-fish measurement, statistics, and the
publication figure (including the representative image panels) directly from the
raw ND2 z-stacks. Point `raw_data/` at the acquisition folder (see the setup cell)
and run all cells.

## Data availability

The raw confocal ND2 z-stacks are not tracked in git (size). They are available
at **[deposit DOI / repository — to be added]** and support full reproduction,
including the representative image panels. The committed per-fish measurement
tables are sufficient to reproduce all reported statistics and the quantification
plot.

## License

Released under the MIT License (see `LICENSE`).
