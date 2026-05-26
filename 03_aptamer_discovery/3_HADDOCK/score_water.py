#!/usr/bin/env python3
"""
score_water.py
Ranks all 40 HADDOCK Phase 1 pairs using water-refined structures.

Logic:
  - Scores from it1/water/file.list (pre-computed by HADDOCK, 1-based index)
  - Clusters from it1/water/analysis/cluster.out (indices = file.list positions)
  - Cluster scored by top4_avg of member scores
  - Pairs ranked by top cluster top4_avg ascending
  - Cross-referenced with haddock_pairs_manifest.csv for metadata
"""

import re
import csv
from pathlib import Path

HADDOCK_DIR = Path("haddock_runs")
MANIFEST   = Path("haddock_pairs_manifest.csv")
OUT_CSV    = Path("haddock_water_scores.csv")

# ── Load manifest for metadata ───────────────────────────────────────────────
manifest = {}
with open(MANIFEST) as f:
    for row in csv.DictReader(f):
        key = f"{row['BestVariant']}__{row['aptamer']}"
        manifest[key] = {
            "class":      row["Class"],
            "phosphosite": row["Phosphosite"],
            "hdock_score": float(row["HDOCKScore"]),
        }

def parse_file_list(path):
    """Returns list of (score, struct_name) in file.list order (index 0 = position 1)."""
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.match(r'"(?:PREVIT:)?(.+?)"\s+\{\s*([-\d.eE+]+)\s*\}', line)
            if m:
                entries.append((float(m.group(2)), m.group(1)))
    return entries  # already sorted best→worst by HADDOCK

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
                members = numbers[1:]  # rest are members (center excluded)
                clusters.append({
                    "id":         cid,
                    "center_idx": center,
                    "members":    [center] + members,  # include center in scoring
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
                    return int(vals[0].strip())  # first value = AIR violations
                except ValueError:
                    return None
    return None

# ── Main loop ────────────────────────────────────────────────────────────────
results = []

for pair_dir in sorted(HADDOCK_DIR.iterdir()):
    if not pair_dir.is_dir():
        continue

    name      = pair_dir.name
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

    entries  = parse_file_list(file_list)   # list[(score, name)], 0-indexed
    clusters = parse_cluster_out(clust_f)

    if not entries or not clusters:
        print(f"PARSE ERROR: {name}")
        continue

    # Score each cluster
    cluster_results = []
    for clust in clusters:
        member_scores = []
        for idx in clust["members"]:
            if 1 <= idx <= len(entries):
                score, sname = entries[idx - 1]  # convert 1-based → 0-based
                member_scores.append((score, sname))

        if not member_scores:
            continue

        member_scores.sort(key=lambda x: x[0])  # ascending (best first)
        top4     = member_scores[:4]
        top4_avg = round(sum(s for s, _ in top4) / len(top4), 4)
        best_score, best_struct = member_scores[0]

        # Center score
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

    # Rank clusters by top4_avg
    cluster_results.sort(key=lambda x: x["top4_avg"])
    top = cluster_results[0]

    # AIR violations for best structure
    air_viol = get_air_violations(water_dir, top["best_struct"])

    # Flag unreliable results
    flags = []
    if top["size"] < 4:
        flags.append("small_cluster")
    if air_viol is not None and air_viol > 0:
        flags.append(f"AIR_viol={air_viol}")

    meta = manifest.get(name, {})

    results.append({
        "pair":         name,
        "class":        meta.get("class", "?"),
        "phosphosite":  meta.get("phosphosite", "?"),
        "hdock_score":  meta.get("hdock_score", "?"),
        "n_clusters":   len(cluster_results),
        "top_clust_id": top["cluster_id"],
        "clust_size":   top["size"],
        "top4_avg":     top["top4_avg"],
        "best_score":   top["best_score"],
        "best_struct":  top["best_struct"],
        "air_viol":     air_viol if air_viol is not None else "?",
        "flags":        ";".join(flags) if flags else "OK",
    })

    print(f"{name:<45} clusters={len(cluster_results)}  "
          f"top_size={top['size']:>3}  "
          f"top4_avg={top['top4_avg']:>10.4f}  "
          f"best={top['best_score']:>10.4f}  "
          f"flags={';'.join(flags) if flags else 'OK'}")

# ── Final ranking ────────────────────────────────────────────────────────────
results.sort(key=lambda x: x["top4_avg"])

print(f"\n{'='*95}")
print(f"PHASE 1 WATER-REFINED CLUSTER RANKING — ALL {len(results)} PAIRS")
print(f"{'='*95}")
print(f"{'Rk':<4} {'Pair':<42} {'Cl':<3} {'Ph':<18} {'NClust':>6} {'Size':>5} "
      f"{'Top4Avg':>10} {'Best':>10}  Flags")
print("-" * 95)

for i, r in enumerate(results, 1):
    marker = " ◄" if i <= 20 else ""
    print(f"{i:<4} {r['pair']:<42} {r['class']:<3} {r['phosphosite']:<18} "
          f"{r['n_clusters']:>6} {r['clust_size']:>5} "
          f"{r['top4_avg']:>10.4f} {r['best_score']:>10.4f}  "
          f"{r['flags']}{marker}")

print(f"\n{'='*95}")
print(f"TOP 20 SHORTLIST (water-refined, cluster-ranked)")
print(f"{'='*95}")
print(f"{'Rk':<4} {'Pair':<42} {'Class':<6} {'Phosphosite':<18} "
      f"{'Top4Avg':>10} {'Size':>5}  Flags")
print("-" * 95)
for i, r in enumerate(results[:20], 1):
    print(f"{i:<4} {r['pair']:<42} {r['class']:<6} {r['phosphosite']:<18} "
          f"{r['top4_avg']:>10.4f} {r['clust_size']:>5}  {r['flags']}")

# ── Write CSV ────────────────────────────────────────────────────────────────
fields = ["rank","pair","class","phosphosite","hdock_score","n_clusters",
          "top_clust_id","clust_size","top4_avg","best_score","best_struct",
          "air_viol","flags"]
with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for i, r in enumerate(results, 1):
        writer.writerow({"rank": i, **r})

print(f"\nSaved: {OUT_CSV}")
