#!/usr/bin/env python3
"""
compute_specificity.py
Combines pTau and hTau water-refined scores and computes
ΔΔG_specificity = Score(hTau) - Score(pTau) for each aptamer.

Inputs:
  - haddock_water_scores.csv        (pTau, all 40)
  - htau_water_scores.csv   (hTau, top 20)

Output:
  - haddock_specificity_scores.csv  ranked by ΔΔG desc (most specific first)
"""

import csv
from pathlib import Path

PTAU_CSV = Path("haddock_water_scores.csv")
HTAU_CSV = Path("htau_water_scores.csv")
OUT_CSV  = Path("haddock_specificity_scores.csv")

def load_ptau_scores(path):
    """Return dict keyed by aptamer id with pTau data."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pair = row["pair"]          # e.g. pTau_2.9__Aptamer_1152
            if "__" not in pair:
                continue
            aptamer = pair.split("__", 1)[1]   # Aptamer_1152
            data[aptamer] = {
                "ptau_pair":    pair,
                "class":        row.get("class", "?"),
                "phosphosite":  row.get("phosphosite", "?"),
                "ptau_top4":    float(row["top4_avg"]),
                "ptau_clust":   int(row["top_clust_id"]),
                "ptau_size":    int(row["clust_size"]),
                "ptau_flags":   row.get("flags", "OK"),
            }
    return data

def load_htau_scores(path):
    """Return dict keyed by aptamer id with hTau data."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pair = row["pair"]          # e.g. hTau_2__Aptamer_1152
            if "__" not in pair:
                continue
            aptamer = pair.split("__", 1)[1]
            data[aptamer] = {
                "htau_pair":   pair,
                "htau_top4":   float(row["top4_avg"]),
                "htau_clust":  int(row["top_clust_id"]),
                "htau_size":   int(row["clust_size"]),
                "htau_flags":  row.get("flags", "OK"),
            }
    return data

ptau = load_ptau_scores(PTAU_CSV)
htau = load_htau_scores(HTAU_CSV)

results = []

for aptamer, p in ptau.items():
    if aptamer not in htau:
        # This pTau pair was not in the hTau run set (e.g. only top 20 have hTau)
        continue
    h = htau[aptamer]

    ptau_score = p["ptau_top4"]
    htau_score = h["htau_top4"]
    ddg = round(htau_score - ptau_score, 4)   # ΔΔG = hTau - pTau

    results.append({
        "aptamer":      aptamer,
        "class":        p["class"],
        "phosphosite":  p["phosphosite"],
        "ptau_pair":    p["ptau_pair"],
        "htau_pair":    h["htau_pair"],
        "ptau_top4":    ptau_score,
        "htau_top4":    htau_score,
        "ddg_specific": ddg,
        "ptau_clust":   p["ptau_clust"],
        "ptau_size":    p["ptau_size"],
        "ptau_flags":   p["ptau_flags"],
        "htau_clust":   h["htau_clust"],
        "htau_size":    h["htau_size"],
        "htau_flags":   h["htau_flags"],
    })

# Rank by specificity: larger ΔΔG = better (hTau worse than pTau)
results.sort(key=lambda x: x["ddg_specific"], reverse=True)

print("="*120)
print("PHASE 1 SPECIFICITY RANKING (ΔΔG = Score(hTau) - Score(pTau))")
print("="*120)
print(f"{'Rk':<4} {'Aptamer':<14} {'Class':<5} {'Phosphosite':<12} "
      f"{'pTau_top4':>10} {'hTau_top4':>10} {'ΔΔG_spec':>10}  "
      f"pFlags / hFlags")
print("-"*120)

for i, r in enumerate(results, 1):
    print(f"{i:<4} {r['aptamer']:<14} {r['class']:<5} {r['phosphosite']:<12} "
          f"{r['ptau_top4']:>10.4f} {r['htau_top4']:>10.4f} {r['ddg_specific']:>10.4f}  "
          f"{r['ptau_flags']} / {r['htau_flags']}")

fields = ["rank","aptamer","class","phosphosite",
          "ptau_pair","htau_pair",
          "ptau_top4","htau_top4","ddg_specific",
          "ptau_clust","ptau_size","ptau_flags",
          "htau_clust","htau_size","htau_flags"]

with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for i, r in enumerate(results, 1):
        w.writerow({"rank": i, **r})

print(f"\nSaved: {OUT_CSV}")