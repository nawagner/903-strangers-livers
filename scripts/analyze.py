#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "pandas",
#   "scikit-learn",
#   "umap-learn",
# ]
# ///
"""
analyze.py — preprocess the ARCHS4 human liver TSV into JSON assets
for the "903 Strangers' Livers" interactive analysis.

Inputs:
  ../demo/RNA_seq_dataset.zip  (human_liver.tsv inside)

Outputs (analysis/data/):
  meta.json                 — matrix shape, sparsity, library size stats
  samples.json              — 903 GSM IDs + per-sample library sizes
  top_genes.json            — top 50 genes by median log-CPM with liver-marker annotation
  normalization_demo.json   — raw vs log-CPM for two samples across canonical liver genes
  pca.json                  — per-sample PC1..PC5 coordinates + variance explained
  loadings.json             — top +/- gene loadings for PC1 and PC2
  matrix_sample.json        — 30 genes x 50 samples raw-count block for the hook viz

Runs via: uv run analysis/scripts/analyze.py
"""
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import umap

# ---- paths ----
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]            # .../Interesting Group Presentation
DEMO = ROOT / "demo"
OUT = ROOT / "analysis" / "data"
OUT.mkdir(parents=True, exist_ok=True)

ZIP_PATH = DEMO / "RNA_seq_dataset.zip"
TSV_NAME = "human_liver.tsv"
TSV_PATH = Path("/tmp/rnaseq_explore") / TSV_NAME

# ---- 1. load ----
if not TSV_PATH.exists():
    TSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(f"Extracting {ZIP_PATH.name} ...")
    with zipfile.ZipFile(ZIP_PATH) as zf:
        zf.extract(TSV_NAME, TSV_PATH.parent)

print(f"Reading {TSV_PATH} ...")
df = pd.read_csv(TSV_PATH, sep="\t", index_col=0)
df.index = df.index.map(str)
df = df.fillna(0)
n_genes, n_samples = df.shape
print(f"  shape: {n_genes} genes x {n_samples} samples")

samples = list(map(str, df.columns))
genes = list(map(str, df.index))

# ---- 2. normalization: log2(CPM + 1) ----
lib_size = df.sum(axis=0)
cpm = df.divide(lib_size, axis=1) * 1e6
logcpm = np.log2(cpm + 1.0)
print(f"  library size range: {int(lib_size.min()):,} -- {int(lib_size.max()):,}")
print(f"  log2(CPM+1) range: {logcpm.values.min():.2f} -- {logcpm.values.max():.2f}")

# ---- 3. top genes by median log-CPM ----
median_logcpm = logcpm.median(axis=1).sort_values(ascending=False)
top50 = median_logcpm.head(50)

LIVER_MARKERS = {
    "ALB":      "Albumin — the most abundant blood protein, made exclusively in hepatocytes.",
    "TF":       "Transferrin — iron transport, a classic hepatocyte secretion.",
    "HP":       "Haptoglobin — scavenges free hemoglobin; acute-phase liver protein.",
    "SERPINA1": "Alpha-1 antitrypsin — the major blood protease inhibitor, liver-made.",
    "APOB":     "Apolipoprotein B — the backbone of LDL, assembled in hepatocytes.",
    "APOA1":    "Apolipoprotein A-I — the main HDL protein, liver-derived.",
    "APOA2":    "Apolipoprotein A-II — HDL component from liver.",
    "FGA":      "Fibrinogen alpha — coagulation factor synthesized in the liver.",
    "FGB":      "Fibrinogen beta — coagulation factor synthesized in the liver.",
    "FGG":      "Fibrinogen gamma — coagulation factor synthesized in the liver.",
    "TTR":      "Transthyretin — thyroid hormone transport, liver-made.",
    "CYP3A4":   "Cytochrome P450 3A4 — the workhorse drug-metabolizing enzyme of the liver.",
    "ORM1":     "Orosomucoid — acute phase protein.",
    "ORM2":     "Orosomucoid 2 — acute phase protein paralog.",
    "AFP":      "Alpha-fetoprotein — fetal liver marker, re-expressed in hepatocellular carcinoma.",
    "GC":       "Vitamin D binding protein — synthesized by hepatocytes.",
    "SAA1":     "Serum amyloid A1 — acute-phase liver protein.",
    "SAA2":     "Serum amyloid A2 — acute-phase liver protein.",
    "C3":       "Complement C3 — liver-synthesized innate immunity protein.",
    "A2M":      "Alpha-2-macroglobulin — liver-synthesized protease inhibitor.",
    "AHSG":     "Fetuin-A — liver-secreted carrier protein.",
    "ITIH1":    "Inter-alpha-trypsin inhibitor heavy chain 1 — liver-synthesized.",
    "ITIH2":    "Inter-alpha-trypsin inhibitor heavy chain 2 — liver-synthesized.",
    "APCS":     "Serum amyloid P — acute phase, made in the liver.",
    "F2":       "Prothrombin — central clotting factor, liver-synthesized.",
    "AMBP":     "Alpha-1-microglobulin/bikunin precursor — liver-synthesized.",
    "HPX":      "Hemopexin — binds free heme in blood, liver-made.",
    "CP":       "Ceruloplasmin — main copper transport protein, liver-made.",
}

top_gene_rows = []
for gene, med in top50.items():
    values = logcpm.loc[gene].values
    top_gene_rows.append({
        "gene": str(gene),
        "median_logcpm": round(float(med), 3),
        "mean_logcpm": round(float(values.mean()), 3),
        "min_logcpm": round(float(values.min()), 3),
        "max_logcpm": round(float(values.max()), 3),
        "q25": round(float(np.percentile(values, 25)), 3),
        "q75": round(float(np.percentile(values, 75)), 3),
        "is_liver_marker": str(gene) in LIVER_MARKERS,
        "annotation": LIVER_MARKERS.get(str(gene)),
    })

with (OUT / "top_genes.json").open("w") as f:
    json.dump({
        "genes": top_gene_rows,
        "markers_in_top50": sum(1 for r in top_gene_rows if r["is_liver_marker"]),
        "total_markers_checked": len(LIVER_MARKERS),
    }, f, indent=2)
print(f"  wrote top_genes.json ({len(top_gene_rows)} rows, "
      f"{sum(1 for r in top_gene_rows if r['is_liver_marker'])} are liver markers)")

# ---- 4. normalization demo ----
lib_sorted = lib_size.sort_values()
small_sample = lib_sorted.index[int(len(lib_sorted) * 0.05)]
large_sample = lib_sorted.index[int(len(lib_sorted) * 0.95)]

demo_gene_candidates = [
    "ALB", "TF", "HP", "SERPINA1", "APOB", "FGB", "FGG",
    "TTR", "CYP3A4", "AFP", "A2M", "C3",
]
canonical_genes = [g for g in demo_gene_candidates if g in df.index]

norm_rows = []
for gene in canonical_genes:
    norm_rows.append({
        "gene": gene,
        "raw_small": int(df.loc[gene, small_sample]),
        "raw_large": int(df.loc[gene, large_sample]),
        "logcpm_small": round(float(logcpm.loc[gene, small_sample]), 3),
        "logcpm_large": round(float(logcpm.loc[gene, large_sample]), 3),
    })

with (OUT / "normalization_demo.json").open("w") as f:
    json.dump({
        "small_sample": {
            "gsm": str(small_sample),
            "library_size": int(lib_size[small_sample]),
        },
        "large_sample": {
            "gsm": str(large_sample),
            "library_size": int(lib_size[large_sample]),
        },
        "library_size_ratio": round(float(lib_size[large_sample] / lib_size[small_sample]), 2),
        "genes": norm_rows,
    }, f, indent=2)
print(f"  wrote normalization_demo.json "
      f"(ratio {round(float(lib_size[large_sample] / lib_size[small_sample]), 1)}x)")

# ---- 5. PCA on top variable genes ----
var_per_gene = logcpm.var(axis=1)
top_var_genes = var_per_gene.sort_values(ascending=False).head(2000).index
X = logcpm.loc[top_var_genes].values.T  # samples x genes
X_mean = X.mean(axis=0, keepdims=True)
X_std = X.std(axis=0, keepdims=True)
X_std[X_std == 0] = 1.0
Xz = (X - X_mean) / X_std

pca = PCA(n_components=10)
pcs = pca.fit_transform(Xz)
var_explained = pca.explained_variance_ratio_
print(f"  PCA variance explained (top 5): "
      f"{[round(float(v) * 100, 1) for v in var_explained[:5]]}%")

pca_rows = []
for i, gsm in enumerate(samples):
    pca_rows.append({
        "gsm": gsm,
        "pc": [round(float(pcs[i, k]), 3) for k in range(5)],
    })

with (OUT / "pca.json").open("w") as f:
    json.dump({
        "samples": pca_rows,
        "variance_explained": [round(float(v), 5) for v in var_explained],
        "n_genes_used": int(len(top_var_genes)),
        "method": "log2(CPM+1) -> top 2000 variable genes -> z-score per gene -> PCA",
    }, f, indent=2)
print(f"  wrote pca.json")

# ---- 5b. UMAP on the same matrix used for PCA ----
# n_neighbors: 15 is the default. Bumping to 50 tells UMAP to consider
#   substantially more global structure when laying out points, which has
#   two effects we want: (1) it pulls satellite clusters closer to the rest
#   of the embedding so the canvas doesn't have empty space between blobs,
#   and (2) it stops UMAP from over-fragmenting the dense central region.
# min_dist: the minimum distance UMAP will allow between points in the
#   embedding. 0.1 is the default; 0.0 produces the tightest possible
#   packing inside clusters, which is what we want for a compact scatter.
# spread: the effective scale of embedded points. Lower values squeeze
#   the whole layout to a smaller absolute range, complementing min_dist=0.0.
print("running UMAP ...")
UMAP_PARAMS = dict(n_neighbors=50, min_dist=0.0, spread=0.6)
reducer = umap.UMAP(
    n_components=2,
    metric="euclidean",
    random_state=42,
    verbose=False,
    **UMAP_PARAMS,
)
umap_coords = reducer.fit_transform(Xz)
umap_rows = []
for i, gsm in enumerate(samples):
    umap_rows.append({
        "gsm": gsm,
        "u": [round(float(umap_coords[i, 0]), 3), round(float(umap_coords[i, 1]), 3)],
    })
with (OUT / "umap.json").open("w") as f:
    json.dump({
        "samples": umap_rows,
        "n_genes_used": int(len(top_var_genes)),
        **UMAP_PARAMS,
        "method": (
            "log2(CPM+1) -> top 2000 variable genes -> z-score per gene -> "
            f"UMAP (n_neighbors={UMAP_PARAMS['n_neighbors']}, "
            f"min_dist={UMAP_PARAMS['min_dist']})"
        ),
    }, f, indent=2)
print(f"  wrote umap.json")

# ---- 6. loadings ----
loadings_matrix = pca.components_  # (10, 2000)

def top_loadings(component_idx: int, k: int = 20):
    loads = loadings_matrix[component_idx]
    order_pos = np.argsort(-loads)[:k]
    order_neg = np.argsort(loads)[:k]
    pos = [{"gene": str(top_var_genes[i]), "loading": round(float(loads[i]), 4)} for i in order_pos]
    neg = [{"gene": str(top_var_genes[i]), "loading": round(float(loads[i]), 4)} for i in order_neg]
    return pos, neg

pc1_pos, pc1_neg = top_loadings(0)
pc2_pos, pc2_neg = top_loadings(1)

with (OUT / "loadings.json").open("w") as f:
    json.dump({
        "PC1": {"positive": pc1_pos, "negative": pc1_neg},
        "PC2": {"positive": pc2_pos, "negative": pc2_neg},
    }, f, indent=2)
print(f"  wrote loadings.json")

# ---- 7. samples list ----
with (OUT / "samples.json").open("w") as f:
    json.dump({
        "samples": samples,
        "library_size": {s: int(lib_size[s]) for s in samples},
    }, f)
print(f"  wrote samples.json")

# ---- 8. matrix sample for hook visualization ----
rng = np.random.default_rng(42)
hook_genes_top = [str(g) for g in top50.index[:15]]
other_pool = [str(g) for g in top_var_genes if str(g) not in hook_genes_top][:200]
hook_genes_other = list(rng.choice(other_pool, size=15, replace=False))
hook_genes = hook_genes_top + hook_genes_other
hook_samples = list(rng.choice(samples, size=50, replace=False))

raw_block = df.loc[hook_genes, hook_samples]
matrix_sample = {
    "genes": list(map(str, hook_genes)),
    "samples": list(map(str, hook_samples)),
    "raw": [[int(x) for x in row] for row in raw_block.values.tolist()],
}
with (OUT / "matrix_sample.json").open("w") as f:
    json.dump(matrix_sample, f)
print(f"  wrote matrix_sample.json ({len(hook_genes)}x{len(hook_samples)})")

# ---- 9. overall metadata ----
zero_frac = float((df.values == 0).sum()) / (n_genes * n_samples)
with (OUT / "meta.json").open("w") as f:
    json.dump({
        "n_genes": int(n_genes),
        "n_samples": int(n_samples),
        "total_cells": int(n_genes * n_samples),
        "raw_max": int(df.values.max()),
        "raw_min": int(df.values.min()),
        "zero_fraction": round(zero_frac, 4),
        "lib_size_min": int(lib_size.min()),
        "lib_size_max": int(lib_size.max()),
        "lib_size_median": int(lib_size.median()),
    }, f, indent=2)
print(f"  wrote meta.json")
print("done.")
