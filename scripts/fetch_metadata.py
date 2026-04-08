#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "requests",
# ]
# ///
"""
fetch_metadata.py — pull GEO metadata for the 903 ARCHS4 liver samples via
NCBI E-utilities (esearch + esummary on db=gds).

Politeness:
  - Batched queries (50 accessions per esearch via OR)
  - Sleeps between requests to stay under 3 req/sec (NCBI default)
  - Caches per-sample to data/meta_cache.json; resumable on rerun

Outputs:
  analysis/data/samples_meta.json — {gsm: {gse, title, source, characteristics, summary}}
  analysis/data/meta_stats.json   — high-level summary used by the site

Run: uv run analysis/scripts/fetch_metadata.py
"""
import json
import re
import sys
import time
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
DATA = ROOT / "analysis" / "data"
CACHE = DATA / "meta_cache.json"
OUT = DATA / "samples_meta.json"
STATS = DATA / "meta_stats.json"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH = 50
SLEEP = 0.4  # seconds between requests

# ---- load sample list from analyze.py output ----
samples_json = DATA / "samples.json"
if not samples_json.exists():
    sys.exit(f"{samples_json} not found. Run analyze.py first.")
with samples_json.open() as f:
    gsm_list = json.load(f)["samples"]
print(f"Loaded {len(gsm_list)} GSM accessions")

# ---- load cache ----
cache: dict = {}
if CACHE.exists():
    try:
        with CACHE.open() as f:
            cache = json.load(f)
        print(f"Resuming from cache: {len(cache)} already fetched")
    except Exception:
        cache = {}

todo = [g for g in gsm_list if g not in cache]
print(f"Still to fetch: {len(todo)}")


def save_cache():
    tmp = CACHE.with_suffix(".tmp")
    with tmp.open("w") as f:
        json.dump(cache, f)
    tmp.replace(CACHE)


def esearch_batch(accessions: list[str]) -> dict[str, str]:
    """Return {GSM accession: UID}."""
    term = " OR ".join(f"{a}[ACCN]" for a in accessions)
    params = {
        "db": "gds",
        "term": term,
        "retmode": "json",
        "retmax": str(len(accessions) * 2),
    }
    r = requests.get(f"{EUTILS}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    uids = r.json().get("esearchresult", {}).get("idlist", [])
    return uids  # we'll map via esummary below


def esummary_batch(uids: list[str]) -> dict[str, dict]:
    if not uids:
        return {}
    params = {
        "db": "gds",
        "id": ",".join(uids),
        "retmode": "json",
    }
    r = requests.get(f"{EUTILS}/esummary.fcgi", params=params, timeout=60)
    r.raise_for_status()
    result = r.json().get("result", {})
    out = {}
    for uid in result.get("uids", []):
        out[uid] = result.get(uid, {})
    return out


def parse_record(rec: dict) -> dict:
    """Keep only the useful fields from an NCBI GDS record."""
    title = rec.get("title", "") or ""
    summary = rec.get("summary", "") or ""
    accession = rec.get("accession", "") or ""
    gpl = rec.get("gpl", "") or ""
    taxon = rec.get("taxon", "") or ""

    # GDS esummary bundles the parent GSE under 'gse' (string, possibly with prefix stripped)
    gse_field = rec.get("gse", "") or ""
    gse = ""
    if gse_field:
        gse = f"GSE{gse_field}" if not str(gse_field).startswith("GSE") else str(gse_field)

    # extract characteristics from 'extrelations' or from the title/summary heuristically
    source_name = rec.get("source_name", "") or ""

    # heuristic disease/tissue keywords — we'll refine in-site
    text = f"{title} {summary} {source_name}".lower()
    disease_tag = "unknown"
    if re.search(r"\bhepatocellular carcinoma\b|\bhcc\b", text):
        disease_tag = "HCC"
    elif re.search(r"\bcholangiocarcinoma\b|\biCC\b", text):
        disease_tag = "cholangiocarcinoma"
    elif re.search(r"\bfibrosis\b|\bcirrhosis\b", text):
        disease_tag = "fibrosis/cirrhosis"
    elif re.search(r"\bnafld\b|\bnash\b|\bsteatohepatitis\b|\bsteatosis\b|\bfatty liver\b", text):
        disease_tag = "NAFLD/NASH"
    elif re.search(r"\bhepatitis\b|\bhcv\b|\bhbv\b", text):
        disease_tag = "hepatitis"
    elif re.search(r"\bhepg2\b|\bhuh7\b|\bhuh-7\b|\bcell line\b", text):
        disease_tag = "cell line"
    elif re.search(r"\bnormal\b|\bhealthy\b|\bcontrol\b", text):
        disease_tag = "normal"

    return {
        "accession": accession,
        "title": title[:300],
        "summary": summary[:600],
        "gse": gse,
        "gpl": gpl,
        "taxon": taxon,
        "source_name": source_name,
        "disease_tag": disease_tag,
    }


def fetch_batch(batch: list[str]):
    # step 1: esearch gives us UIDs for this batch
    uids = esearch_batch(batch)
    time.sleep(SLEEP)
    if not uids:
        # if the batch failed completely, mark each as empty so we don't retry forever
        for g in batch:
            cache[g] = {"accession": g, "error": "esearch_empty"}
        return
    # step 2: esummary on those UIDs
    summaries = esummary_batch(uids)
    time.sleep(SLEEP)

    # map summaries back to GSM accessions
    by_acc = {}
    for uid, rec in summaries.items():
        parsed = parse_record(rec)
        acc = parsed["accession"]
        if acc:
            by_acc[acc] = parsed

    for g in batch:
        if g in by_acc:
            cache[g] = by_acc[g]
        else:
            cache[g] = {"accession": g, "error": "not_in_summary"}


# ---- main loop ----
t0 = time.time()
for i in range(0, len(todo), BATCH):
    batch = todo[i:i + BATCH]
    try:
        fetch_batch(batch)
    except Exception as e:
        print(f"  batch {i}: error {e} — will retry up to 2x")
        ok = False
        for attempt in range(2):
            time.sleep(2 + attempt * 3)
            try:
                fetch_batch(batch)
                ok = True
                break
            except Exception as e2:
                print(f"    retry {attempt + 1}: {e2}")
        if not ok:
            for g in batch:
                cache[g] = {"accession": g, "error": "fetch_failed"}
    save_cache()
    done = sum(1 for g in gsm_list if g in cache)
    elapsed = time.time() - t0
    rate = (i + len(batch)) / max(elapsed, 0.001)
    eta = (len(todo) - (i + len(batch))) / max(rate, 0.001)
    print(f"  {done}/{len(gsm_list)} fetched  ({rate:.1f}/s, ETA {eta:.0f}s)")

# ---- write final per-sample and aggregate stats ----
final = {g: cache.get(g, {"accession": g, "error": "missing"}) for g in gsm_list}
with OUT.open("w") as f:
    json.dump(final, f)
print(f"wrote {OUT.name}")

# aggregate stats
disease_counts: dict[str, int] = {}
gse_counts: dict[str, int] = {}
for rec in final.values():
    d = rec.get("disease_tag", "unknown")
    disease_counts[d] = disease_counts.get(d, 0) + 1
    gse = rec.get("gse") or "unknown"
    gse_counts[gse] = gse_counts.get(gse, 0) + 1

with STATS.open("w") as f:
    json.dump({
        "n_samples_fetched": sum(1 for r in final.values() if not r.get("error")),
        "n_samples_failed": sum(1 for r in final.values() if r.get("error")),
        "disease_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
        "gse_counts": dict(sorted(gse_counts.items(), key=lambda kv: -kv[1])),
        "n_gse": len([k for k in gse_counts if k != "unknown"]),
    }, f, indent=2)
print(f"wrote {STATS.name}")
print(f"done in {time.time() - t0:.0f}s")
