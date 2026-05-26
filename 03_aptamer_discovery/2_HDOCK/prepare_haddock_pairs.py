#!/usr/bin/env python3
"""
prepare_haddock_pairs.py
Selects best pair per aptamer strictly from IsHit=True confirmed pairs.
This guarantees AIR files exist for every pair selected.
"""
import json
import pandas as pd
from pathlib import Path

SHORTLIST_XLSX   = Path("wave1_shortlist_outputs_20260318_085331/wave1_shortlist.shortlist.xlsx")
HADDOCK_INPUTS   = Path("wave1_shortlist_outputs_20260318_085331/haddock_inputs/")
OUTPUT_DIR       = Path("./haddock_pairs/")
PHOSPHO_JSON     = Path("./phospho_sites.json")

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    # Load sheets
    df_selected  = pd.read_excel(SHORTLIST_XLSX, sheet_name='Selected_Aptamers')
    df_hpairs    = pd.read_excel(SHORTLIST_XLSX, sheet_name='HADDOCK_Pairs')  # IsHit=True only
    df_jobs      = pd.read_excel(SHORTLIST_XLSX, sheet_name='PerJob_Scores')

    with open(PHOSPHO_JSON) as f:
        phospho_map = json.load(f)

    # For each aptamer, pick the best variant STRICTLY from confirmed hit pairs
    # Sorted by HdockScore ascending (most negative = best)
    df_hpairs_sorted = df_hpairs.sort_values('HdockScore', ascending=True)
    best_hit_per_apt = df_hpairs_sorted.drop_duplicates(subset='aptamer', keep='first')
    best_hit_map     = best_hit_per_apt.set_index('aptamer')

    print(f"{'='*70}")
    print(f"  Phase 1 — 40 Best Confirmed-Hit Pairs for HADDOCK")
    print(f"{'='*70}")
    print(f"  {'#':<4} {'Aptamer':<16} {'Cls':<4} {'BestHitVariant':<16} "
          f"{'Phosphosite':<16} {'HDOCKScore':<12} {'AIR'}")
    print(f"  {'-'*70}")

    pairs = []
    tier_map = dict(zip(df_selected['aptamer'], df_selected['Tier']))
    cls_map  = dict(zip(df_selected['aptamer'], df_selected['Class']))

    for i, row in df_selected.iterrows():
        aptamer = row['aptamer']

        if aptamer not in best_hit_map.index:
            print(f"  {i+1:<4} {aptamer:<16} *** No confirmed hit found — skipping ***")
            continue

        best_hit   = best_hit_map.loc[aptamer]
        best_var   = best_hit['variant']
        best_score = best_hit['HdockScore']

        phospho_res = phospho_map.get(best_var, [])
        phospho_str = '/'.join([
            f"pT{r}" if r in [181, 212, 217, 231] else f"pS{r}"
            for r in phospho_res
        ])

        # Get best pose file from PerJob_Scores
        match = df_jobs[
            (df_jobs['aptamer'] == aptamer) &
            (df_jobs['variant'] == best_var)
        ]
        best_pose = match['BestPoseFile'].values[0] if not match.empty else "N/A"

        # Verify AIR files exist and are non-empty
        base      = f"{best_var}__{aptamer}"
        air_files = {
            'rec_actives':  HADDOCK_INPUTS / f"{base}__rec_actives.txt",
            'rec_passives': HADDOCK_INPUTS / f"{base}__rec_passives.txt",
            'apt_actives':  HADDOCK_INPUTS / f"{base}__apt_actives.txt",
            'apt_passives': HADDOCK_INPUTS / f"{base}__apt_passives.txt",
        }
        air_ok  = all(p.exists() and p.stat().st_size > 0 for p in air_files.values())
        air_str = "✅" if air_ok else "❌"

        pairs.append({
            'rank':         len(pairs) + 1,
            'aptamer':      aptamer,
            'Class':        cls_map.get(aptamer, '?'),
            'BestVariant':  best_var,
            'Phosphosite':  phospho_str,
            'HDOCKScore':   best_score,
            'Tier':         tier_map.get(aptamer, '?'),
            'BestPoseFile': best_pose,
            'AIR_rec_act':  str(air_files['rec_actives']),
            'AIR_rec_pas':  str(air_files['rec_passives']),
            'AIR_apt_act':  str(air_files['apt_actives']),
            'AIR_apt_pas':  str(air_files['apt_passives']),
            'AIR_OK':       air_ok,
        })

        print(f"  {len(pairs):<4} {aptamer:<16} {cls_map.get(aptamer,'?'):<4} "
              f"{best_var:<16} {phospho_str:<16} {best_score:<12} {air_str}")

    # Save manifest
    df_out = pd.DataFrame(pairs)
    manifest_path = OUTPUT_DIR / "haddock_pairs_manifest.csv"
    df_out.to_csv(manifest_path, index=False)

    print(f"\n{'='*70}")
    print(f"  Total pairs          : {len(pairs)}")
    print(f"  AIR files confirmed  : {df_out['AIR_OK'].sum()} / {len(pairs)}")

    print(f"\n  Phosphosite coverage:")
    for ps, count in df_out['Phosphosite'].value_counts().items():
        print(f"    {ps:<22} : {count} aptamers")

    print(f"\n  Design class distribution:")
    for cls, count in df_out['Class'].value_counts().sort_index().items():
        print(f"    Class {cls:<16} : {count} aptamers")

    if not df_out['AIR_OK'].all():
        print(f"\n  ❌ Still missing AIR files for:")
        for _, r in df_out[~df_out['AIR_OK']].iterrows():
            print(f"    {r['BestVariant']}__{r['aptamer']}")
    else:
        print(f"\n  ✅ All 40 AIR file sets confirmed! Ready for HADDOCK.")

    print(f"\n  Manifest: {manifest_path}")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
