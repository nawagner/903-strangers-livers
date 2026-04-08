# 903 Strangers' Livers

An interactive exploratory analysis of the ARCHS4 human-liver RNA-seq dataset
([demo/RNA_seq_dataset.zip](../demo/RNA_seq_dataset.zip)), built as a pudding-cool-style
scrollytelling web piece for the ASAP webinar *"Best Practices for AI in Scientific Research."*

The point of the piece is not to publish a finding — it's to show what happens
when you hand a pile of unlabeled biological data to a researcher (or an AI agent)
and let it ask the data to explain itself.

## What's here

```
analysis/
├── index.html          # the scrollytelling piece
├── styles.css          # ASAP brand palette, Noto Sans + Noto Serif, layout
├── main.js             # D3 charts, scroll observer, PCA brushing, recolor
├── assets/
│   └── liver.svg       # public-domain liver illustration (Häggström, Wikimedia Commons)
├── data/               # pre-computed JSON (committed; regenerate via scripts/)
│   ├── meta.json
│   ├── samples.json
│   ├── top_genes.json
│   ├── normalization_demo.json
│   ├── pca.json
│   ├── umap.json
│   ├── loadings.json
│   ├── matrix_sample.json
│   ├── samples_meta.json        # per-sample GEO metadata, enriched
│   ├── meta_stats.json          # disease / study breakdown
│   ├── meta_cache.json          # fetch_metadata.py cache (resumable)
│   └── study_cache.json         # enrich_metadata.py cache (resumable)
└── scripts/
    ├── analyze.py               # TSV → JSON (normalization, PCA, gene ranking)
    ├── fetch_metadata.py        # NCBI E-utilities sample-level metadata
    └── enrich_metadata.py       # add GSE study-level context + disease tagging
```

The browser never touches the 79 MB TSV. All heavy compute runs offline in
Python and the web client just reads the small JSON files under `data/`.

## Running locally

```bash
# (from the analysis/ directory)
python3 -m http.server 8765
open http://localhost:8765/
```

Append `?visible` to the URL to bypass the scroll-triggered reveal animations —
useful for printing, screenshot capture, or previewing the full piece at a glance.

## Regenerating the data

All three Python scripts run via [uv](https://github.com/astral-sh/uv) with
inline dependencies declared in a PEP 723 header — no venv setup required.

```bash
# (from the project root, one directory above analysis/)
uv run analysis/scripts/analyze.py           # ~30s — TSV → JSON
uv run analysis/scripts/fetch_metadata.py    # ~30s — 903 GSMs via NCBI E-utilities
uv run analysis/scripts/enrich_metadata.py   # ~2min — 60 GSE studies + re-tagging
```

Each fetcher caches its results to disk and resumes on rerun. The NCBI calls are
rate-limited to stay under the default 3 requests/second limit.

## How the analysis works

**Normalization.** Each sample's raw counts are divided by its library size,
multiplied by a million, and log-transformed (`log2(CPM + 1)`). The library
sizes in this dataset span 106k to 132M reads — a 1,236× range — so this step
is not optional.

**Top-genes ranking.** Genes are ranked by median `log2(CPM+1)` across all 903
samples. The top 50 are annotated against a hand-curated list of canonical
hepatocyte-secretome markers (ALB, TF, HP, SERPINA1, APOB, fibrinogens, etc.).
Eleven of the top 50 turn out to be liver markers — strong corroboration that
the filename means what it says.

**PCA.** The 2,000 most variable genes (by variance of log-CPM) are z-scored
and fed to scikit-learn's PCA. PC1 captures 36.0% of total variance; PC2
captures 13.7% — so a single flat scatter captures half the structure.

**UMAP.** The same z-scored matrix is also passed through `umap-learn` with
`n_neighbors=50, min_dist=0.0, spread=0.6, random_state=42`, producing a
non-linear 2D projection that preserves local neighbourhoods at the cost
of global distances. The high `n_neighbors` value pulls satellite clusters
closer to the rest of the embedding so the canvas isn't dominated by
whitespace; `min_dist=0.0` and `spread=0.6` produce a tight, compact layout.
The piece renders both projections so the geometric contrast is visible.

**Metadata fetching.** [fetch_metadata.py](scripts/fetch_metadata.py) pulls
per-sample GEO records via `esearch` + `esummary` on `db=gds`. That gives us
each sample's GSE (study of origin). [enrich_metadata.py](scripts/enrich_metadata.py)
then pulls the *study-level* record for each of the 60 unique GSE accessions
— those records carry the rich summary text.

**Sample tagging.** Tagging happened in three escalating passes:

1. **Regex pass** (`enrich_metadata.py`): pattern-matching on per-sample
   text. It was wrong on 53 of 60 studies — for example, a 7-sample "normal"
   study that turned out to be tumor-educated platelets.
2. **LLM-per-study pass** (`apply_llm_tags.py`): three parallel Claude Code
   sub-agents read every study and applied a controlled vocabulary. Better,
   but propagating one tag per study is wrong when a study has multiple
   sample groups (e.g. an iPSC liver-bud paper that also contributes its
   primary fetal-hepatocyte and adult-hepatocyte references — those
   reference samples were getting tagged "iPSC/ESC differentiation").
3. **LLM-per-sample-group pass** (`apply_group_tags.py`, current):
   [build_sample_group_audit.py](scripts/build_sample_group_audit.py) bins
   the 903 samples into **86 unique groups** by `(gse, sample_summary,
   sample_source_name)`. Three parallel sub-agents classify each group with
   per-sample evidence taking precedence over study context. Result:
   `data/llm_group_labels.json`. The diff against the previous study-level
   pass is in [data/group_tag_diff.json](data/group_tag_diff.json) (57 of 86
   groups changed).

The final hand-curated breakdown:

| Tag                                  | Samples |
| ------------------------------------ | ------- |
| normal fetal liver                   | 217     |
| iPSC/ESC-derived hepatocyte-like     | 178     |
| primary hepatocytes (in vitro)       |  99     |
| non-tumor adjacent liver             |  79     |
| normal adult liver                   |  77     |
| HCC                                  |  52     |
| viral hepatitis                      |  46     |
| hepatic organoid                     |  30     |
| NAFLD/NASH                           |  26     |
| hepatic stellate cell                |  17     |
| cholangiocarcinoma                   |  14     |
| liver-derived hematopoietic          |  14     |
| non-liver                            |  13     |
| cholestasis                          |  11     |
| hepatoblastoma                       |  10     |
| liver metastasis (non-liver primary) |   9     |
| hepatocyte cell line                 |   7     |
| liver-resident immune cell           |   4     |

Note how the picture changes between passes:

- The big single-study mistake (GSE96981 with 299 samples) was previously
  lumped into "iPSC/ESC differentiation". The per-group pass correctly
  splits its 299 samples into 206 "Human_fetal_liver" and 93
  "Human_adult_liver" — they were never iPSC samples, they were the
  primary-hepatocyte references for the iPSC paper.
- "non-liver" dropped from 42 to 13 samples once we stopped propagating
  study tags to multi-condition studies.
- A new category, **non-tumor adjacent liver** (79 samples), captures the
  matched-control samples from cancer studies — biologically these are
  near-normal liver, but they were previously getting tagged as the cancer
  type of the parent study.

The single biggest takeaway is still the same: this is **not** a clean
slice of adult human liver. The largest single category is fetal liver
tissue (24% of the dataset), and another 31% are some kind of in-vitro
model system (hepatocyte cultures, organoids, iPSC-derived cells, cell
lines). The "903 human liver samples" framing is honest only if you know
exactly which 903 you got.

## Chart inventory

1. **Hero grid** — 903 animated dots drifting into place; each dot is one sample.
2. **Raw-matrix heatmap** — 30 genes × 50 samples, log-scale coloring, hover tooltips.
3. **Top 50 genes** — horizontal bar chart with IQR whiskers, liver markers highlighted in emerald.
4. **Normalization slider** — dumbbell plot across 12 canonical liver genes, transitioning between raw counts → CPM → `log2(CPM+1)`.
5. **Act I PCA** — 903-sample scatter, drag-to-brush with a side panel showing the selection's centroid and aligned loading genes.
6. **Loadings panel** — top ±15 gene loadings for PC1 and PC2.
7. **Metadata breakdown** — Act II; inferred disease tag counts + per-study sample counts.
8. **Act II PCA** — recolorable scatter with three modes: *no colour / by study (GSE) / by disease tag*. This is the reveal that shows PCA clusters map more tightly to study-of-origin than to biology — i.e. batch effects dominate.
9. **UMAP comparison** — the same 903 samples × 2,000 genes projected non-linearly with UMAP, with the same three recolor modes. Cluster boundaries are sharper but inter-cluster distances become meaningless — the standard PCA-vs-UMAP trade-off.

## Citation

The processed counts come from ARCHS4. If you use this analysis or build on it,
please cite the ARCHS4 paper:

> Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC,
> Ma'ayan A. **Massive mining of publicly available RNA-seq data from human and
> mouse.** *Nature Communications* 9, 1366 (2018).
> [doi.org/10.1038/s41467-018-03751-6](https://doi.org/10.1038/s41467-018-03751-6)

The 903-sample human-liver subset itself was redistributed by Alexander Lachmann
on [Kaggle](https://www.kaggle.com/datasets/lachmann12/human-liver-rnaseq-gene-expression-903-samples)
under CC0. Sample-level metadata comes from
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) via E-utilities, fetched live by
[fetch_metadata.py](scripts/fetch_metadata.py).

## Brand

Palette and typography are drawn from [parkinsonsroadmap.org](https://parkinsonsroadmap.org/)
and the presentation's existing [brand reference](../presentation/CLAUDE.md).
Fonts are Noto Sans and Noto Serif, served from Google Fonts. The liver mark in
the header is by Mikael Häggström, M.D., released to the public domain via
[Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Liver.svg).

## Dependencies

- Python 3.10+, `uv`
- `numpy`, `pandas`, `scikit-learn`, `umap-learn` (declared inline in [analyze.py](scripts/analyze.py))
- `requests` (declared inline in [fetch_metadata.py](scripts/fetch_metadata.py) and [enrich_metadata.py](scripts/enrich_metadata.py))
- D3 v7 and Observable Plot, loaded from jsdelivr ESM — **no build step**
- A modern browser with ES modules and `IntersectionObserver` (any of the last 5 years of Chrome/Firefox/Safari)

## Deploying to Vercel

This is a fully static site with no build step. From the `analysis/` directory:

```bash
vercel --prod
```

Vercel auto-detects it as a static project and serves the files directly.
