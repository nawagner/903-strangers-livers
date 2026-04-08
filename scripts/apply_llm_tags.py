#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
apply_llm_tags.py — apply LLM-generated study labels to samples_meta.json
and rewrite meta_stats.json so the visualization picks them up.

Inputs:
  data/llm_study_labels.json   — array of {gse, primary_tag, secondary_tag,
                                            is_liver_focused, confidence, rationale}
                                  produced by aggregating sub-agent outputs
  data/samples_meta.json       — current per-sample records (will be rewritten)

Outputs:
  data/samples_meta.json       — disease_tag replaced with LLM primary_tag,
                                  plus new fields: llm_primary_tag,
                                  llm_secondary_tag, llm_confidence,
                                  llm_is_liver_focused, llm_rationale
  data/meta_stats.json         — recomputed from new tags
  data/tag_diff.json           — per-study record of regex_tag vs llm_tag,
                                  for the comparison report
"""
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analysis" / "data"

LLM = json.loads((DATA / "llm_study_labels.json").read_text())
samples = json.loads((DATA / "samples_meta.json").read_text())

# index by gse
by_gse = {entry["gse"]: entry for entry in LLM}

# regex tag for diff report (preserved from current samples_meta before we overwrite)
old_regex_tags: dict[str, str] = {}
for gsm, rec in samples.items():
    gse_field = (rec.get("gse") or "").split(";")[0].strip() or "(no_gse)"
    old_regex_tags.setdefault(gse_field, rec.get("disease_tag", "unknown"))

# apply
applied = 0
unmatched = []
for gsm, rec in samples.items():
    gse_field = (rec.get("gse") or "").split(";")[0].strip()
    if not gse_field:
        unmatched.append(gsm)
        continue
    label = by_gse.get(gse_field)
    if not label:
        unmatched.append(gsm)
        continue
    rec["disease_tag"] = label["primary_tag"]
    rec["llm_primary_tag"] = label["primary_tag"]
    rec["llm_secondary_tag"] = label.get("secondary_tag")
    rec["llm_is_liver_focused"] = label.get("is_liver_focused", True)
    rec["llm_confidence"] = label.get("confidence", "medium")
    rec["llm_rationale"] = label.get("rationale", "")
    samples[gsm] = rec
    applied += 1

(DATA / "samples_meta.json").write_text(json.dumps(samples))
print(f"applied LLM tags to {applied}/{len(samples)} samples ({len(unmatched)} unmatched)")

# rebuild meta_stats.json from the new tags
disease_counts: Counter = Counter()
gse_counts: Counter = Counter()
liver_focused_counts: Counter = Counter()
for rec in samples.values():
    disease_counts[rec.get("disease_tag", "unknown")] += 1
    gse_counts[rec.get("gse") or "unknown"] += 1
    liver_focused_counts[bool(rec.get("llm_is_liver_focused", True))] += 1

stats = {
    "n_samples_fetched": sum(1 for r in samples.values() if not r.get("error")),
    "n_samples_failed": sum(1 for r in samples.values() if r.get("error")),
    "disease_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
    "gse_counts": dict(sorted(gse_counts.items(), key=lambda kv: -kv[1])),
    "n_gse": len([k for k in gse_counts if k != "unknown"]),
    "tagging_method": "LLM (GPT-class sub-agents reading per-GSE study text)",
    "n_liver_focused": liver_focused_counts.get(True, 0),
    "n_non_liver": liver_focused_counts.get(False, 0),
}
(DATA / "meta_stats.json").write_text(json.dumps(stats, indent=2))
print(f"rewrote meta_stats.json")
print(f"  liver-focused samples: {stats['n_liver_focused']}/{len(samples)}")
print(f"  non-liver samples:     {stats['n_non_liver']}/{len(samples)}")

# diff report (per study, showing where regex disagreed with LLM)
diff_rows = []
for entry in LLM:
    gse = entry["gse"]
    old = old_regex_tags.get(gse, "(no regex tag)")
    new = entry["primary_tag"]
    n = sum(1 for r in samples.values() if (r.get("gse") or "").split(";")[0].strip() == gse)
    diff_rows.append({
        "gse": gse,
        "n_samples": n,
        "regex_tag": old,
        "llm_tag": new,
        "agreed": old == new,
        "is_liver_focused": entry.get("is_liver_focused", True),
        "confidence": entry.get("confidence", "medium"),
        "rationale": entry.get("rationale", ""),
    })

(DATA / "tag_diff.json").write_text(json.dumps(diff_rows, indent=2))
print(f"wrote tag_diff.json")

# inline summary of disagreements
disagreements = [r for r in diff_rows if not r["agreed"]]
print()
print(f"=== {len(disagreements)} of {len(diff_rows)} studies relabeled ===")
print(f"{'GSE':<14} {'n':>4}  {'regex':<32} {'->':<3} {'LLM':<32}")
print("-" * 90)
for r in sorted(disagreements, key=lambda r: -r["n_samples"]):
    print(f"{r['gse']:<14} {r['n_samples']:>4}  {r['regex_tag']:<32} -> {r['llm_tag']:<32}")

print()
print("=== final disease tag distribution (LLM) ===")
for tag, count in sorted(disease_counts.items(), key=lambda kv: -kv[1]):
    print(f"  {tag:<30} {count:>4}")
