#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
build_sample_group_audit.py — produce a per-group audit table.

The previous tagging assumed every sample in a study shared the study's
biology. That's wrong: ARCHS4 keyword-searches at the sample level, so
multi-condition studies can contribute heterogeneous samples to the
"liver" subset (e.g. a study about iPSC-derived liver buds also
contributes the primary fetal/adult hepatocyte reference samples).

This script groups samples by (gse, sample_summary, source_name) and
produces a table of unique groups for per-group LLM reclassification.

Outputs:
  data/sample_group_audit.json   one row per unique (gse, summary, source)
                                  with example sample titles and current tag
"""
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analysis" / "data"

samples = json.loads((DATA / "samples_meta.json").read_text())
study_cache = json.loads((DATA / "study_cache.json").read_text())

groups: dict = defaultdict(list)
for gsm, r in samples.items():
    gse = (r.get("gse") or "").split(";")[0].strip() or "(no_gse)"
    summary = (r.get("summary") or "").strip()
    source = (r.get("source_name") or "").strip()
    groups[(gse, summary, source)].append({
        "gsm": gsm,
        "title": (r.get("title") or "").strip(),
    })

rows = []
for (gse, summary, source), members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    study = study_cache.get(gse, {})
    sample_titles = sorted({m["title"] for m in members if m["title"]})
    titles_str = " | ".join(sample_titles[:6])
    if len(sample_titles) > 6:
        titles_str += f" | (+{len(sample_titles) - 6} more)"

    # representative sample for current_tag (since all members of a group
    # share the same study, this is uniform within a group)
    first = members[0]["gsm"]
    current_tag = samples[first].get("disease_tag", "unknown")

    rows.append({
        "group_id": f"{gse}__{len(rows):03d}",
        "gse": gse,
        "n_samples": len(members),
        "current_tag": current_tag,
        "sample_summary": summary or "(empty)",
        "sample_source_name": source or "(empty)",
        "sample_titles_sample": titles_str,
        "study_title": (study.get("title") or "").strip(),
        "study_summary": (study.get("summary") or "").strip()[:500],
        "gsm_list": [m["gsm"] for m in members],
    })

(DATA / "sample_group_audit.json").write_text(json.dumps(rows, indent=2))
print(f"wrote sample_group_audit.json — {len(rows)} unique groups, "
      f"{sum(r['n_samples'] for r in rows)} samples")

# inline summary of largest groups
print()
print(f"{'group':<18} {'gse':<14} {'n':>4} {'current_tag':<28} {'sample_summary'}")
print("-" * 130)
for r in rows[:30]:
    s = r["sample_summary"][:50]
    print(f"{r['group_id']:<18} {r['gse']:<14} {r['n_samples']:>4} {r['current_tag']:<28} {s}")
print(f"... +{max(0, len(rows) - 30)} more")
