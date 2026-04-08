#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
build_study_audit.py — produce a per-GSE audit table for the 60 studies
that contributed samples to the human-liver subset.

Reads:
  data/samples_meta.json
  data/study_cache.json

Writes:
  data/study_audit.csv      one row per unique GSE
  data/study_audit.json     same content, structured for LLM consumption

Each row contains everything needed to re-classify a study by reading:
  gse, n_samples, current_regex_tag, title, summary, overall_design,
  taxon, sample_titles_sample, all_disease_tags_seen.
"""
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analysis" / "data"

samples_meta = json.loads((DATA / "samples_meta.json").read_text())
study_cache = json.loads((DATA / "study_cache.json").read_text())

# Group samples by their primary GSE (first one if there's a semicolon list)
samples_by_gse: dict[str, list[dict]] = defaultdict(list)
for gsm, rec in samples_meta.items():
    gse_field = (rec.get("gse") or "").strip()
    if not gse_field:
        gse_field = "(no_gse)"
    primary = re.split(r"[;,]", gse_field)[0].strip() or "(no_gse)"
    samples_by_gse[primary].append({"gsm": gsm, **rec})

rows = []
for gse, sample_list in sorted(samples_by_gse.items(), key=lambda kv: -len(kv[1])):
    study = study_cache.get(gse, {})
    n = len(sample_list)
    # current regex tag(s) seen — should be uniform within a study, but show counts in case
    tag_counter = Counter(s.get("disease_tag", "unknown") for s in sample_list)
    primary_tag = tag_counter.most_common(1)[0][0] if tag_counter else "unknown"

    # gather a few sample titles to give the reviewer a feel for what's in the study
    sample_titles = sorted(set(s.get("title", "") for s in sample_list if s.get("title")))
    sample_titles_str = " | ".join(sample_titles[:8])
    if len(sample_titles) > 8:
        sample_titles_str += f" | (+{len(sample_titles) - 8} more)"

    rows.append({
        "gse": gse,
        "n_samples": n,
        "current_regex_tag": primary_tag,
        "study_title": (study.get("title") or "").strip(),
        "study_summary": (study.get("summary") or "").strip(),
        "overall_design": (study.get("overall_design") or "").strip(),
        "taxon": (study.get("taxon") or "").strip(),
        "sample_titles_sample": sample_titles_str,
        "all_tags_in_study": dict(tag_counter),
    })

# CSV — flatten all_tags_in_study to a string
csv_path = DATA / "study_audit.csv"
with csv_path.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "gse", "n_samples", "current_regex_tag",
        "study_title", "study_summary", "overall_design",
        "taxon", "sample_titles_sample", "all_tags_in_study",
    ])
    for r in rows:
        writer.writerow([
            r["gse"],
            r["n_samples"],
            r["current_regex_tag"],
            r["study_title"],
            r["study_summary"],
            r["overall_design"],
            r["taxon"],
            r["sample_titles_sample"],
            ";".join(f"{k}={v}" for k, v in r["all_tags_in_study"].items()),
        ])
print(f"wrote {csv_path.name}  ({len(rows)} rows)")

# JSON — for LLM ingestion
json_path = DATA / "study_audit.json"
json_path.write_text(json.dumps(rows, indent=2))
print(f"wrote {json_path.name}")

# Inline summary
print()
print(f"{'GSE':<20} {'n':>4}  {'regex_tag':<28} {'title (truncated)'}")
print("-" * 110)
for r in rows:
    title = r["study_title"][:60]
    print(f"{r['gse']:<20} {r['n_samples']:>4}  {r['current_regex_tag']:<28} {title}")
