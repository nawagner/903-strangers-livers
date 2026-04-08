#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "requests",
# ]
# ///
"""
enrich_metadata.py — enrich samples_meta.json with study-level context.

The per-sample esummary records we already fetched return minimal fields
(often just "Liver_2"). The richer disease/tissue context lives at the
*study* (GSE) level. This script:

  1. loads samples_meta.json
  2. collects unique GSE accessions
  3. queries NCBI E-utilities for study-level summaries
  4. re-tags each sample by running keyword matching on the study text
  5. rewrites samples_meta.json and meta_stats.json
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

SAMPLES_META = DATA / "samples_meta.json"
STATS_OUT = DATA / "meta_stats.json"
STUDY_CACHE = DATA / "study_cache.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SLEEP = 0.4


def classify(text: str) -> str:
    t = text.lower()
    if re.search(r"\bhepatocellular carcinoma\b|\bhcc\b|\bliver cancer\b|\bliver tumou?r\b", t):
        return "HCC"
    if re.search(r"\bcholangiocarcinoma\b|\bintrahepatic cholangio\b|\bbile duct cancer\b", t):
        return "cholangiocarcinoma"
    if re.search(r"\bfibrosis\b|\bcirrhosis\b|\bcirrhotic\b", t):
        return "fibrosis/cirrhosis"
    if re.search(r"\bnafld\b|\bnash\b|\bsteatohepatitis\b|\bsteatosis\b|\bfatty liver\b", t):
        return "NAFLD/NASH"
    if re.search(r"\bhepatitis\b|\bhcv\b|\bhbv\b|\bviral liver\b", t):
        return "hepatitis"
    if re.search(r"\bhepg2\b|\bhuh7\b|\bhuh-7\b|\bhepatoma cell\b|\bcell line\b|\bhep3b\b", t):
        return "cell line"
    if re.search(r"\bhepatic stellate\b|\bstellate cell\b", t):
        return "stellate cell"
    if re.search(r"\bhepatocyte\b", t) and re.search(r"\bprimary\b|\bisolated\b", t):
        return "primary hepatocytes"
    if re.search(r"\borganoid\b", t):
        return "organoid"
    if re.search(r"\bembryonic\b|\bfetal\b|\bipsc\b|\bpluripotent\b|\bdifferentiation\b|\bstem cell\b", t):
        return "stem cell / differentiation"
    if re.search(r"\bnormal\b|\bhealthy\b|\bcontrol\b", t):
        return "normal"
    return "unknown"


def fetch_study(gse: str) -> dict:
    """Fetch a GSE study record via esearch+esummary."""
    # esearch
    r = requests.get(
        f"{EUTILS}/esearch.fcgi",
        params={"db": "gds", "term": f"{gse}[ACCN]", "retmode": "json"},
        timeout=30,
    )
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"gse": gse, "title": "", "summary": "", "error": "not_found"}
    time.sleep(SLEEP)

    # esummary
    r = requests.get(
        f"{EUTILS}/esummary.fcgi",
        params={"db": "gds", "id": ",".join(ids), "retmode": "json"},
        timeout=60,
    )
    r.raise_for_status()
    result = r.json().get("result", {})
    # pick the record whose accession starts with "GSE" (not GPL/GDS/etc.)
    best = None
    for uid in result.get("uids", []):
        rec = result.get(uid, {})
        acc = rec.get("accession", "")
        if acc == gse:
            best = rec
            break
    if best is None:
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            if rec.get("accession", "").startswith("GSE"):
                best = rec
                break
    if best is None:
        return {"gse": gse, "title": "", "summary": "", "error": "no_gse_record"}
    return {
        "gse": gse,
        "title": best.get("title", "") or "",
        "summary": best.get("summary", "") or "",
        "overall_design": best.get("overalldesign", "") or "",
        "taxon": best.get("taxon", "") or "",
        "n_samples": best.get("n_samples", 0),
    }


def main():
    if not SAMPLES_META.exists():
        sys.exit(f"{SAMPLES_META} not found. Run fetch_metadata.py first.")

    with SAMPLES_META.open() as f:
        samples = json.load(f)
    print(f"loaded {len(samples)} samples")

    # collect unique GSEs; samples may have composite accessions like "GSE80126;81080"
    all_gse = set()
    for rec in samples.values():
        gse = rec.get("gse") or ""
        if not gse:
            continue
        for g in re.split(r"[;,]", gse):
            g = g.strip()
            if g.startswith("GSE"):
                all_gse.add(g)
    print(f"unique studies to fetch: {len(all_gse)}")

    # resumable cache
    studies: dict = {}
    if STUDY_CACHE.exists():
        try:
            studies = json.loads(STUDY_CACHE.read_text())
        except Exception:
            studies = {}

    todo = sorted(g for g in all_gse if g not in studies)
    print(f"todo: {len(todo)}")
    t0 = time.time()
    for i, gse in enumerate(todo):
        try:
            studies[gse] = fetch_study(gse)
        except Exception as e:
            print(f"  {gse}: error {e}")
            studies[gse] = {"gse": gse, "title": "", "summary": "", "error": str(e)}
        time.sleep(SLEEP)
        if (i + 1) % 10 == 0 or i == len(todo) - 1:
            STUDY_CACHE.write_text(json.dumps(studies, indent=2))
            print(f"  {i + 1}/{len(todo)} ({(i + 1) / max(time.time() - t0, 0.001):.1f}/s)")
    STUDY_CACHE.write_text(json.dumps(studies, indent=2))

    # re-tag samples using study text + original sample text
    disease_counts: dict[str, int] = {}
    gse_counts: dict[str, int] = {}
    for gsm, rec in samples.items():
        sample_gse = (rec.get("gse") or "").split(";")[0].strip()
        study = studies.get(sample_gse, {})
        text = " ".join([
            rec.get("title") or "",
            rec.get("summary") or "",
            rec.get("source_name") or "",
            study.get("title") or "",
            study.get("summary") or "",
            study.get("overall_design") or "",
        ])
        tag = classify(text)
        rec["disease_tag"] = tag
        rec["study_title"] = study.get("title", "")
        rec["study_summary"] = (study.get("summary") or "")[:400]
        samples[gsm] = rec
        disease_counts[tag] = disease_counts.get(tag, 0) + 1
        gse = rec.get("gse") or "unknown"
        gse_counts[gse] = gse_counts.get(gse, 0) + 1

    with SAMPLES_META.open("w") as f:
        json.dump(samples, f)
    print(f"rewrote {SAMPLES_META.name}")

    stats = {
        "n_samples_fetched": sum(1 for r in samples.values() if not r.get("error")),
        "n_samples_failed": sum(1 for r in samples.values() if r.get("error")),
        "disease_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
        "gse_counts": dict(sorted(gse_counts.items(), key=lambda kv: -kv[1])),
        "n_gse": len([k for k in gse_counts if k != "unknown"]),
    }
    with STATS_OUT.open("w") as f:
        json.dump(stats, f, indent=2)
    print(f"rewrote {STATS_OUT.name}")
    print("disease breakdown:")
    for tag, n in sorted(disease_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {tag:30s} {n}")


if __name__ == "__main__":
    main()
