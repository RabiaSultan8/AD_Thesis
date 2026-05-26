#!/usr/bin/env python3
"""
score_htau_water.py
Post-HADDOCK analysis for healthy Tau (hTau) dockings.

- Mirrors score_water.py used for Phase 1 pTau runs.
- Parses it1/water/file.list and analysis/cluster.out for each hTau run.
- Ranks pairs by top4_avg of best cluster (ascending).
- Writes htau_water_scores.csv.
"""

import re
import csv
from pathlib import Path

HTAU_DIR = Path("htau_runs")                    # root with all hTau_*__Aptamer_* dirs
OUT_CSV  = Path("htau_water_scores.csv")

def parse_file_list(path):
    """Returns list of (score, struct_name) in file.list order (index 0 = position 1)."""
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.match(r'"(?:PREVIT:)?(.+?)"\s+\{\s*([-+\d.eE]+)\s*\}', line)
            if m:
                entries.append((float(m.group(2)), m.group(1)))
    return entries

def parse_cluster_out(path):
    """
    Returns list of clusters sorted by cluster ID:
      [{"id": N, "center_idx": int, "members": [int, ...]}, ...]
    Indices are 1-based positions in file.list.
    """
    clusters = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            m = re.match(r'Cluster\s+(\d+)\s+->\s+([\d\s]+)', line)
            if m:
                cid     = int(m.group(1))
                numbers = [int(x) for x in m.group(2).split()]
                center  = numbers[0]
                members = numbers[1:]
                clusters.append({
                    "id":         cid,
                    "center_idx": center,
                    "members":    [center] + members,
                })
    return clusters

def get_air_violations(water_dir, struct_name):
    """Read AIR violations from the best structure's REMARK."""
    pdb = water_dir / struct_name
    if not pdb.exists():
        return None
    with open(pdb) as f:
        for line in f:
            if line.startswith("REMARK violations.:"):
                vals = line.split(":")[1].strip().split(",")
                try:
                    return int(vals[0].strip())
                except ValueError:
                    return None
    return None

results = []

for pair_dir in sorted(HTAU_DIR.iterdir()):
    if not pair_dir.is_dir():
        continue

    name      = pair_dir.name          # e.g. hTau_2__Aptamer_1152
    water_dir = pair_dir / "run1" / "structures" / "it1" / "water"
    ana_dir   = water_dir / "analysis"
    file_list = water_dir / "file.list"
    clust_f   = ana_dir   / "cluster.out"

    if not file_list.exists():
        print(f"MISSING file.list: {name}")
        continue
    if not clust_f.exists():
        print(f"MISSING cluster.out: {name}")
        continue

    entries  = parse_file_list(file_list)
    clusters = parse_cluster_out(clust_f)

    if not entries or not clusters:
        print(f"PARSE ERROR: {name}")
        continue

    cluster_results = []
    for clust in clusters:
        member_scores = []
        for idx in clust["members"]:
            if 1 <= idx <= len(entries):
                score, sname = entries[idx - 1]
                member_scores.append((score, sname))
        if not member_scores:
            continue

        member_scores.sort(key=lambda x: x[0])  # ascending
        top4     = member_scores[:4]
        top4_avg = round(sum(s for s, _ in top4) / len(top4), 4)
        best_score, best_struct = member_scores[0]

        ci = clust["center_idx"]
        center_score = entries[ci - 1][0] if 1 <= ci <= len(entries) else None

        cluster_results.append({
            "cluster_id":   clust["id"],
            "size":         len(member_scores),
            "top4_avg":     top4_avg,
            "best_score":   round(best_score, 4),
            "best_struct":  best_struct,
            "center_score": round(center_score, 4) if center_score else None,
        })

    if not cluster_results:
        print(f"NO CLUSTERS SCORED: {name}")
        continue

    cluster_results.sort(key=lambda x: x["top4_avg"])
    top = cluster_results[0]

    air_viol = get_air_violations(water_dir, top["best_struct"])

    flags = []
    if top["size"] < 4:
        flags.append("small_cluster")
    if air_viol is not None and air_viol > 0:
        flags.append(f"AIR_viol={air_viol}")

    results.append({
        "pair":         name,
        "n_clusters":   len(cluster_results),
        "top_clust_id": top["cluster_id"],
        "clust_size":   top["size"],
        "top4_avg":     top["top4_avg"],
        "best_score":   top["best_score"],
        "best_struct":  top["best_struct"],
        "air_viol":     air_viol if air_viol is not None else "?",
        "flags":        ";".join(flags) if flags else "OK",
    })

    print(f"{name:<30} clusters={len(cluster_results):>3}  "
          f"top_size={top['size']:>3}  "
          f"top4_avg={top['top4_avg']:>10.4f}  "
          f"best={top['best_score']:>10.4f}  "
          f"flags={';'.join(flags) if flags else 'OK'}")

# rank ascending by top4_avg
results.sort(key=lambda x: x["top4_avg"])

print("\n" + "="*95)
print(f"HADDOCK Phase 1 — hTau water-refined cluster ranking ({len(results)} pairs)")
print("="*95)
print(f"{'Rk':<4} {'Pair':<32} {'NClust':>6} {'Size':>5} "
      f"{'Top4Avg':>10} {'Best':>10}  Flags")
print("-"*95)

for i, r in enumerate(results, 1):
    print(f"{i:<4} {r['pair']:<32} {r['n_clusters']:>6} {r['clust_size']:>5} "
          f"{r['top4_avg']:>10.4f} {r['best_score']:>10.4f}  {r['flags']}")

fields = ["rank","pair","n_clusters","top_clust_id","clust_size",
          "top4_avg","best_score","best_struct","air_viol","flags"]
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for i, r in enumerate(results, 1):
        w.writerow({"rank": i, **r})

print(f"\nSaved: {OUT_CSV}")