#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
apply_group_tags.py — apply per-group LLM labels to samples_meta.json.

Where the previous apply_llm_tags.py propagated one tag per *study*, this
script propagates per-(gse, sample_summary, sample_source_name) groups,
which is the right granularity for ARCHS4 keyword searches.

Inputs:
  data/sample_group_audit.json   the per-group audit table (with gsm lists)
  data/llm_group_labels.json     LLM tag per group_id

Outputs:
  data/samples_meta.json         disease_tag rewritten per-group
  data/meta_stats.json           recomputed
  data/group_tag_diff.json       per-group diff vs prior tag
"""
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analysis" / "data"

groups = json.loads((DATA / "sample_group_audit.json").read_text())
labels = json.loads((DATA / "llm_group_labels.json").read_text())
samples = json.loads((DATA / "samples_meta.json").read_text())

label_by_id = {l["group_id"]: l for l in labels}

# Sanity check that we have a label for every group
missing = [g["group_id"] for g in groups if g["group_id"] not in label_by_id]
if missing:
    raise SystemExit(f"missing labels for {len(missing)} groups: {missing[:5]}...")

applied = 0
diff_rows = []
for g in groups:
    label = label_by_id[g["group_id"]]
    new_tag = label["primary_tag"]
    old_tag = g["current_tag"]
    for gsm in g["gsm_list"]:
        if gsm not in samples:
            continue
        rec = samples[gsm]
        rec["disease_tag"] = new_tag
        rec["llm_primary_tag"] = new_tag
        rec["llm_is_liver_focused"] = label.get("is_liver_focused", True)
        rec["llm_confidence"] = label.get("confidence", "medium")
        rec["llm_rationale"] = label.get("rationale", "")
        rec["llm_group_id"] = g["group_id"]
        samples[gsm] = rec
        applied += 1
    diff_rows.append({
        "group_id": g["group_id"],
        "gse": g["gse"],
        "n_samples": g["n_samples"],
        "old_tag": old_tag,
        "new_tag": new_tag,
        "agreed": old_tag == new_tag,
        "is_liver_focused": label.get("is_liver_focused", True),
        "confidence": label.get("confidence", "medium"),
        "rationale": label.get("rationale", ""),
    })

(DATA / "samples_meta.json").write_text(json.dumps(samples))
print(f"applied per-group tags to {applied}/{len(samples)} samples")

# rebuild stats
disease_counts: Counter = Counter()
gse_counts: Counter = Counter()
liver_focused = Counter()
for rec in samples.values():
    disease_counts[rec.get("disease_tag", "unknown")] += 1
    gse_counts[rec.get("gse") or "unknown"] += 1
    liver_focused[bool(rec.get("llm_is_liver_focused", True))] += 1

stats = {
    "n_samples_fetched": sum(1 for r in samples.values() if not r.get("error")),
    "n_samples_failed": sum(1 for r in samples.values() if r.get("error")),
    "disease_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
    "gse_counts": dict(sorted(gse_counts.items(), key=lambda kv: -kv[1])),
    "n_gse": len([k for k in gse_counts if k != "unknown"]),
    "tagging_method": "LLM (per-(gse,sample_summary) groups via Claude Code sub-agents)",
    "n_liver_focused": liver_focused.get(True, 0),
    "n_non_liver": liver_focused.get(False, 0),
    "n_groups": len(groups),
}
(DATA / "meta_stats.json").write_text(json.dumps(stats, indent=2))
print(f"rewrote meta_stats.json")
print(f"  liver-focused samples: {stats['n_liver_focused']}/{len(samples)}")
print(f"  non-liver samples:     {stats['n_non_liver']}/{len(samples)}")

(DATA / "group_tag_diff.json").write_text(json.dumps(diff_rows, indent=2))
print(f"wrote group_tag_diff.json")

# inline diff report
disagreements = [r for r in diff_rows if not r["agreed"]]
print()
print(f"=== {len(disagreements)} of {len(diff_rows)} groups relabeled ===")
print(f"{'group':<18} {'n':>4}  {'old':<32} {'->':<3} {'new':<35}")
print("-" * 100)
for r in sorted(disagreements, key=lambda r: -r["n_samples"]):
    old = r["old_tag"][:30]
    new = r["new_tag"][:35]
    print(f"{r['group_id']:<18} {r['n_samples']:>4}  {old:<32} -> {new:<35}")

print()
print("=== final disease tag distribution ===")
for tag, count in sorted(disease_counts.items(), key=lambda kv: -kv[1]):
    print(f"  {tag:<40} {count:>4}")
