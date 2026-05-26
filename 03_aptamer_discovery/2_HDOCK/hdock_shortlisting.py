#!/usr/bin/env python3
"""
hdock_shortlisting.py
======================
Post-docking aptamer shortlisting pipeline for HDOCK blind docking results.
Selects top aptamer candidates for HADDOCK flexible docking using a
two-tier, aptamer-centric global ranking strategy.

SCIENTIFIC BASIS FOR THRESHOLDS
---------------------------------
Direct contact cutoff  : 4.5 A
    Huang & Zou (2014) Nucleic Acids Res 42:e55
    HDOCK's own RNA scoring function definition

Proximity cutoff       : 6.0 A
    Dominguez et al. (2003) J Am Chem Soc 125:1731
    HADDOCK AIR (Ambiguous Interaction Restraint) interface definition

SELECTION STRATEGY
-------------------
Tier 1 — Broad binders (diagnostically general)
    Aptamers hitting >= MIN_HITS_TIER1 variants at proximity cutoff
    Ranked by mean HDOCK score across hit variants
    Biologically: recognizes phospho-epitope robustly across pTau landscape

Tier 2 — Specific binders (diagnostically targeted)
    Aptamers hitting < MIN_HITS_TIER1 variants but with exceptional
    HDOCK score (top TOP_PERCENT_TIER2 % of all scores observed)
    Biologically: high-affinity aptamer for specific phospho-site
    e.g. pT217 — most sensitive early AD biomarker
    (Palmqvist et al. 2020 JAMA; Barthelemy et al. 2020 Nat Med)

CSS FORMULA
------------
CSS = 0.40*H_norm + 0.25*E_norm + 0.15*C + 0.20*D_norm
    H_norm : normalized HDOCK score     (primary, validated scoring fn)
    E_norm : normalized PoseScore       (RNA chemistry quality)
    C      : consensus fraction         (pose reproducibility)
    D_norm : normalized proximity       (distance to phospho atoms)

POSE SCORE (E)
---------------
PoseScore = 3.0*BaseContacts + 1.0*SugarContacts
          - 1.5*BackboneContacts + 1.0*GuanineContacts

RNA ATOM CLASSIFICATION
------------------------
Backbone : P, OP1/O1P, OP2/O2P, OP3/O3P, O5', O3'
Sugar    : C1'-C5', O4', O2' (ribose ring)
Base     : all remaining nucleobase atoms
Guanine  : base atoms where resname in {G, DG, GUA, DGU, RG}

HARD REJECTION CRITERIA
-------------------------
no_proximity      : BestMinDist > 6.0 A across all poses
too_few_contacts  : NumContacts < 2 in best pose
backbone_dominant : BackboneContacts > 60% of total in best pose

INPUTS
-------
    docking_results_wave1/
    pTau_1.4/
        Aptamer_47_dock.out
        Aptamer_47/
            model_1.pdb  ...  model_10.pdb
    phospho_sites.json  ->  {"pTau_1.4": [217], ...}

OUTPUTS
--------
    {prefix}_outputs_{timestamp}/
    shortlist.xlsx
        GlobalRanking     <- all aptamers ranked globally
        Tier1_Broad       <- broad binders
        Tier2_Specific    <- specific binders
        Selected_Aptamers <- combined final selection
        HADDOCK_Pairs     <- (aptamer, variant) pairs for HADDOCK
        PerJob_Scores     <- per (aptamer, variant) CSS scores
        HitMatrix         <- 257x16 proximity hit data
    all_poses.csv
    hit_matrix.csv        <- for figure generation
    haddock_inputs/
        pTau_1.4__Aptamer_47__rec_actives.txt
        ...

USAGE
------
    python hdock_shortlisting.py \\
        --input_dir docking_results_wave1 \\
        --phospho_json phospho_sites.json \\
        --output_prefix wave1_shortlist \\
        --top_n_poses 10 \\
        --processes 32
"""

import os
import re
import json
import logging
import argparse
import numpy as np
import pandas as pd
from math import inf
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# ══════════════════════════════════════════════════════════════
#   PARAMETERS
# ══════════════════════════════════════════════════════════════

DIRECT_CUTOFF  = 4.5    # A — Huang & Zou 2014 NAR
PROX_CUTOFF    = 6.0    # A — Dominguez et al. 2003 JACS / HADDOCK AIR

HARD_MIN_CONTACTS      = 2
HARD_MAX_BACKBONE_FRAC = 0.60

MIN_HITS_TIER1    = 3     # variants hit to qualify as broad binder
TOP_PERCENT_TIER2 = 10    # top X% HDOCK score for specific binder
MAX_HADDOCK_INPUT = 40    # hard cap on aptamers going to HADDOCK

W_HDOCK   = 0.40
W_ENERGY  = 0.25
W_CONSENS = 0.15
W_DIST    = 0.20

W_BASE     =  3.0
W_SUGAR    =  1.0
W_BACKBONE = -1.5
W_GUANINE  =  1.0

PHOS_ATOMS     = {"P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3"}
BACKBONE_ATOMS = PHOS_ATOMS | {"O5'", "O3'", "O5*", "O3*"}
SUGAR_ATOMS    = {
    "C1'", "C2'", "C3'", "C4'", "C5'", "O4'",
    "C1*", "C2*", "C3*", "C4*", "C5*", "O4*", "O2'", "O2*",
}
GUANINE_RESNAMES = {"G", "DG", "GUA", "DGU", "RG"}


# ══════════════════════════════════════════════════════════════
#   HDOCK .out PARSER
# ══════════════════════════════════════════════════════════════

def parse_hdock_out(out_file):
    """
    Parses HDOCK .out and returns (hdock_score, confidence) for top pose.
    HDOCK .out columns: tx ty tz rx ry rz Score Lscore Confidence
    Score is negative — more negative = stronger binding.
    Returns (None, None) on failure.
    """
    try:
        with open(out_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if re.match(r'^[A-Za-z/]', line):
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        score = float(parts[6])
                        conf  = float(parts[8]) if len(parts) >= 9 else None
                        return score, conf
                    except (ValueError, IndexError):
                        continue
    except Exception as e:
        logging.debug(f"Could not parse {out_file}: {e}")
    return None, None


# ══════════════════════════════════════════════════════════════
#   PDB PARSING
# ══════════════════════════════════════════════════════════════

def split_complex_pdb(pdb_text):
    """
    Splits createpl complex PDB into receptor and ligand ATOM lines.
    Three-stage fallback: HEADER markers -> TER record -> chain ID.
    """
    lines = pdb_text.splitlines()

    # Strategy 1: HEADER markers
    if any(re.search(r'HEADER\s+(rec|lig)', l, re.IGNORECASE) for l in lines):
        rec, lig, mode = [], [], 'rec'
        for ln in lines:
            if re.search(r'HEADER\s+lig', ln, re.IGNORECASE):
                mode = 'lig'; continue
            if re.search(r'HEADER\s+rec', ln, re.IGNORECASE):
                mode = 'rec'; continue
            if ln.startswith(('ATOM', 'HETATM')):
                (rec if mode == 'rec' else lig).append(ln)
        if rec and lig:
            return rec, lig

    # Strategy 2: TER-based
    rec, lig, past_ter = [], [], False
    for ln in lines:
        if ln.startswith('TER'):
            past_ter = True; continue
        if ln.startswith(('ATOM', 'HETATM')):
            (lig if past_ter else rec).append(ln)
    if rec and lig:
        return rec, lig

    # Strategy 3: Chain ID
    rec, lig = [], []
    for ln in lines:
        if ln.startswith(('ATOM', 'HETATM')):
            chain = ln[21:22].strip() if len(ln) > 21 else 'A'
            (rec if chain in ('', 'A') else lig).append(ln)
    return rec, lig


def parse_atom_line(line):
    """
    Parses PDB ATOM/HETATM line with fixed-width + token fallback.
    Returns dict or None.
    """
    try:
        return dict(
            serial=int(line[6:11]), name=line[12:16].strip(),
            resname=line[17:20].strip(), resid=int(line[22:26]),
            xyz=np.array([float(line[30:38]),
                          float(line[38:46]),
                          float(line[46:54])]),
        )
    except Exception:
        pass
    try:
        t = line.split()
        return dict(serial=int(t[1]), name=t[2], resname=t[3],
                    resid=int(t[4]),
                    xyz=np.array([float(t[5]), float(t[6]), float(t[7])]))
    except Exception:
        return None


# ══════════════════════════════════════════════════════════════
#   RNA ATOM CLASSIFICATION
# ══════════════════════════════════════════════════════════════

def classify_rna_atom(name, resname):
    """
    Classifies RNA atom into backbone/sugar/base/guanine.
    Returns (is_backbone, is_sugar, is_base, is_guanine).
    """
    is_bb  = name in BACKBONE_ATOMS
    is_sug = (not is_bb) and (name in SUGAR_ATOMS)
    is_bas = (not is_bb) and (not is_sug)
    is_gua = is_bas and (resname in GUANINE_RESNAMES
                         or resname.upper().startswith('G'))
    return is_bb, is_sug, is_bas, is_gua


# ══════════════════════════════════════════════════════════════
#   POSE METRICS
# ══════════════════════════════════════════════════════════════

def compute_pose_metrics(rec_lines, lig_lines, phospho_residues):
    """
    Computes geometry and chemistry metrics for one docked pose.
    Returns dict with MinDist, NumContacts, contact type counts,
    PoseScore, and contacts list.
    """
    phos_atoms, phos_xyz = [], []
    for ln in rec_lines:
        a = parse_atom_line(ln)
        if a and a['resid'] in phospho_residues and a['name'] in PHOS_ATOMS:
            phos_atoms.append(a)
            phos_xyz.append(a['xyz'])

    empty = dict(MinDist=inf, NumContacts=0, BaseContacts=0,
                 SugarContacts=0, BackboneContacts=0,
                 GuanineContacts=0, PoseScore=0.0, contacts=[])
    if not phos_xyz:
        return empty

    phos_xyz = np.vstack(phos_xyz)

    apt_atoms, apt_xyz = [], []
    for ln in lig_lines:
        a = parse_atom_line(ln)
        if a:
            apt_atoms.append(a)
            apt_xyz.append(a['xyz'])
    if not apt_xyz:
        return empty

    apt_xyz  = np.vstack(apt_xyz)
    diff     = apt_xyz[:, None, :] - phos_xyz[None, :, :]
    dists    = np.linalg.norm(diff, axis=2)
    min_dist = float(dists.min())

    # Deduplicated contacts: one entry per aptamer atom (closest phos)
    seen = {}
    for ai, pi in np.argwhere(dists <= DIRECT_CUTOFF):
        d = float(dists[ai, pi])
        if ai not in seen or d < seen[ai]['distance']:
            aa, pa = apt_atoms[ai], phos_atoms[pi]
            seen[ai] = dict(
                apt_serial=int(aa['serial']), apt_name=aa['name'],
                apt_resname=aa['resname'],    apt_resid=int(aa['resid']),
                phos_name=pa['name'],         phos_resname=pa['resname'],
                phos_resid=int(pa['resid']),  distance=d,
            )
    contacts = list(seen.values())

    base_c = sugar_c = backbone_c = guanine_c = 0
    for c in contacts:
        is_bb, is_sug, is_bas, is_gua = classify_rna_atom(
            c['apt_name'], c['apt_resname'])
        if is_bb:    backbone_c += 1
        elif is_sug: sugar_c    += 1
        elif is_bas:
            base_c += 1
            if is_gua: guanine_c += 1

    pose_score = (W_BASE * base_c + W_SUGAR * sugar_c
                + W_BACKBONE * backbone_c + W_GUANINE * guanine_c)

    return dict(MinDist=min_dist, NumContacts=len(contacts),
                BaseContacts=base_c, SugarContacts=sugar_c,
                BackboneContacts=backbone_c, GuanineContacts=guanine_c,
                PoseScore=pose_score, contacts=contacts)


# ══════════════════════════════════════════════════════════════
#   WORKER — PER (VARIANT, APTAMER) PAIR
# ══════════════════════════════════════════════════════════════

def analyze_aptamer(args):
    """Multiprocessing worker. Analyzes all poses for one pair."""
    variant, aptamer_dir, out_file, phospho_residues, top_n_poses = args
    aptamer = aptamer_dir.name

    hdock_score, hdock_conf = None, None
    if out_file and Path(out_file).exists():
        hdock_score, hdock_conf = parse_hdock_out(Path(out_file))

    model_files = sorted(aptamer_dir.glob("model_*.pdb"))
    if not model_files:
        model_files = sorted(aptamer_dir.glob("*.pdb"))
    model_files = model_files[:top_n_poses]

    if not model_files:
        return [dict(variant=variant, aptamer=aptamer, pose_index=-1,
                     hdock_score=hdock_score, hdock_conf=hdock_conf,
                     error='no_model_files')]

    records = []
    for mf in model_files:
        digits   = re.findall(r'\d+', mf.stem)
        pose_idx = int(digits[-1]) if digits else -1
        try:
            txt = mf.read_text()
        except Exception as e:
            records.append(dict(variant=variant, aptamer=aptamer,
                                pose_index=pose_idx, error=f'read:{e}',
                                hdock_score=hdock_score,
                                hdock_conf=hdock_conf))
            continue

        rec_lines, lig_lines = split_complex_pdb(txt)
        if not rec_lines or not lig_lines:
            records.append(dict(variant=variant, aptamer=aptamer,
                                pose_index=pose_idx, error='split_failed',
                                hdock_score=hdock_score,
                                hdock_conf=hdock_conf))
            continue

        m = compute_pose_metrics(rec_lines, lig_lines, phospho_residues)
        records.append(dict(
            variant=variant, aptamer=aptamer,
            pose_file=str(mf), pose_index=pose_idx,
            hdock_score=hdock_score, hdock_conf=hdock_conf,
            MinDist=m['MinDist'], NumContacts=m['NumContacts'],
            BaseContacts=m['BaseContacts'], SugarContacts=m['SugarContacts'],
            BackboneContacts=m['BackboneContacts'],
            GuanineContacts=m['GuanineContacts'],
            PoseScore=m['PoseScore'],
            contacts=json.dumps(m['contacts']),
            error=None,
        ))
    return records


# ══════════════════════════════════════════════════════════════
#   PER-JOB AGGREGATION
# ══════════════════════════════════════════════════════════════

def normalize(series, invert=False):
    """Min-max normalize Series to [0,1]. Invert if needed."""
    s  = series.copy().astype(float).replace([inf, -inf], np.nan)
    mn, mx = s.min(), s.max()
    if mx - mn < 1e-8:
        return s.fillna(0.5).clip(0, 1)
    n = (s - mn) / (mx - mn)
    return (1 - n if invert else n).fillna(0.0)


def aggregate_per_job(df_poses):
    """One summary row per (variant, aptamer) pair."""
    rows = []
    for (variant, aptamer), g in df_poses.groupby(['variant', 'aptamer']):
        n_poses   = len(g)
        n_prox    = int((g['MinDist'] <= PROX_CUTOFF).sum())
        cons_frac = n_prox / n_poses if n_poses > 0 else 0.0

        g_s   = g.sort_values(['PoseScore', 'MinDist'],
                               ascending=[False, True]).reset_index(drop=True)
        best  = g_s.iloc[0]
        bd    = float(g['MinDist'].min())
        ms    = float(g['PoseScore'].max())
        bnc   = int(best['NumContacts'])
        bbc   = int(best['BackboneContacts'])
        bbf   = bbc / bnc if bnc > 0 else 0.0

        reject, reason = False, None
        if bd > PROX_CUTOFF:
            reject, reason = True, 'no_proximity'
        elif bnc < HARD_MIN_CONTACTS:
            reject, reason = True, 'too_few_contacts'
        elif bbf > HARD_MAX_BACKBONE_FRAC:
            reject, reason = True, 'backbone_dominant'

        rows.append(dict(
            variant=variant, aptamer=aptamer,
            HdockScore=best.get('hdock_score'),
            HdockConf=best.get('hdock_conf'),
            NumPoses=n_poses, NumPosesProx=n_prox,
            ConsensusFraction=round(cons_frac, 4),
            BestMinDist=round(bd, 3),
            MaxPoseScore=round(ms, 4),
            BestNumContacts=bnc,
            BestBaseContacts=int(best['BaseContacts']),
            BestSugarContacts=int(best['SugarContacts']),
            BestBackboneContacts=bbc,
            BestGuanineContacts=int(best['GuanineContacts']),
            BackboneFraction=round(bbf, 4),
            BestPoseFile=best.get('pose_file', ''),
            BestPoseIndex=int(best.get('pose_index', -1)),
            IsHit=not reject,
            IsDirectHit=(bd <= DIRECT_CUTOFF and not reject),
            HardReject=reject, RejectReason=reason,
        ))
    return pd.DataFrame(rows)


def compute_css(df_jobs):
    """CSS per job, normalized within each variant."""
    df = df_jobs.copy()
    df['CSS'] = np.nan
    for _, g in df.groupby('variant'):
        idx = g.index
        df.loc[idx, 'H_norm'] = normalize(g['HdockScore'], invert=True)
        df.loc[idx, 'E_norm'] = normalize(g['MaxPoseScore'])
        df.loc[idx, 'D_norm'] = normalize(g['BestMinDist'],  invert=True)
        df.loc[idx, 'C']      = g['ConsensusFraction']
        df.loc[idx, 'CSS']    = (
            W_HDOCK * df.loc[idx, 'H_norm']
          + W_ENERGY  * df.loc[idx, 'E_norm']
          + W_CONSENS * df.loc[idx, 'C']
          + W_DIST    * df.loc[idx, 'D_norm']
        ).round(4)
    df.loc[df['HardReject'], 'CSS'] = 0.0
    return df


# ══════════════════════════════════════════════════════════════
#   CROSS-VARIANT AGGREGATION (aptamer-centric)
# ══════════════════════════════════════════════════════════════

def aggregate_per_aptamer(df_jobs):
    """One global summary row per aptamer across all 16 variants."""
    rows = []
    for aptamer, g in df_jobs.groupby('aptamer'):
        hits   = g[g['IsHit']]
        n_hits = len(hits)
        n_dhit = int(g['IsDirectHit'].sum())

        all_sc  = g['HdockScore'].dropna().tolist()
        hit_sc  = hits['HdockScore'].dropna().tolist()

        best_hd  = min(all_sc) if all_sc else None
        mean_hd  = float(np.mean(hit_sc)) if hit_sc else None
        best_css = float(g['CSS'].max())

        hit_variants = sorted(hits['variant'].tolist())
        best_variant = (g.loc[g['HdockScore'].idxmin(), 'variant']
                        if all_sc else None)

        n = int(aptamer.split('_')[1]) if '_' in aptamer else 0
        cls = ('A' if n <= 1250 else 'B' if n <= 2500
               else 'C' if n <= 3750 else 'D')

        rows.append(dict(
            aptamer=aptamer, Class=cls,
            NumHits=n_hits, NumDirectHits=n_dhit,
            NumVariantsTested=len(g),
            MeanHdockHit=round(mean_hd, 3) if mean_hd is not None else None,
            BestHdockScore=round(best_hd, 3) if best_hd is not None else None,
            BestCSS=round(best_css, 4),
            HitVariants='; '.join(hit_variants),
            BestVariant=best_variant,
        ))
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════
#   TWO-TIER SELECTION
# ══════════════════════════════════════════════════════════════

def select_tiers(df_apt, df_jobs):
    """
    Tier 1 — NumHits >= MIN_HITS_TIER1, ranked by MeanHdockHit
    Tier 2 — NumHits < MIN_HITS_TIER1 AND BestHdockScore in top
              TOP_PERCENT_TIER2 % of all scores
    Combined list capped at MAX_HADDOCK_INPUT.
    Returns (tier1, tier2, combined, haddock_pairs).
    """
    df_hits = df_apt[df_apt['NumHits'] >= 1].copy()

    tier1 = (df_hits[df_hits['NumHits'] >= MIN_HITS_TIER1]
             .sort_values('MeanHdockHit', ascending=True)
             .reset_index(drop=True))

    all_scores   = df_hits['BestHdockScore'].dropna()
    score_thresh = float(np.percentile(all_scores, TOP_PERCENT_TIER2))

    tier2 = (df_hits[
        (df_hits['NumHits'] < MIN_HITS_TIER1) &
        (df_hits['BestHdockScore'] <= score_thresh)
    ].sort_values('BestHdockScore', ascending=True)
     .reset_index(drop=True))

    tier1['Tier'] = 'Tier1_Broad'
    tier2['Tier'] = 'Tier2_Specific'

    combined = (pd.concat([tier1, tier2], ignore_index=True)
                  .drop_duplicates(subset='aptamer')
                  .head(MAX_HADDOCK_INPUT))

    # HADDOCK pairs: selected aptamers × their hit variants only
    sel_apts = set(combined['aptamer'])
    haddock_pairs = []
    tier_map = dict(zip(combined['aptamer'], combined['Tier']))
    for _, row in df_jobs[
            df_jobs['aptamer'].isin(sel_apts) &
            df_jobs['IsHit']].iterrows():
        haddock_pairs.append(dict(
            aptamer=row['aptamer'],
            variant=row['variant'],
            Tier=tier_map.get(row['aptamer'], ''),
            HdockScore=row['HdockScore'],
            CSS=row['CSS'],
            BestMinDist=row['BestMinDist'],
            BestPoseFile=row['BestPoseFile'],
        ))

    logging.info(f"Tier 1 (>= {MIN_HITS_TIER1} hits)       : {len(tier1)}")
    logging.info(f"Tier 2 (top {TOP_PERCENT_TIER2}% score) : {len(tier2)}")
    logging.info(f"Combined (cap {MAX_HADDOCK_INPUT})       : {len(combined)}")
    logging.info(f"HADDOCK pairs                            : "
                 f"{len(haddock_pairs)}")

    return tier1, tier2, combined, haddock_pairs


# ══════════════════════════════════════════════════════════════
#   HIT MATRIX
# ══════════════════════════════════════════════════════════════

def build_hit_matrix(df_jobs):
    """
    257x16 hit matrix.
    Values: 0=miss, 1=proximal hit (<=6A), 2=direct hit (<=4.5A)
    """
    apts = sorted(df_jobs['aptamer'].unique(),
                  key=lambda x: int(x.split('_')[1])
                  if '_' in x and x.split('_')[1].isdigit() else 0)
    vars = sorted(df_jobs['variant'].unique())
    mat  = pd.DataFrame(0, index=apts, columns=vars)
    for _, row in df_jobs.iterrows():
        a, v = row['aptamer'], row['variant']
        if a in mat.index and v in mat.columns:
            if row['IsDirectHit']:  mat.loc[a, v] = 2
            elif row['IsHit']:      mat.loc[a, v] = 1
    return mat


# ══════════════════════════════════════════════════════════════
#   OUTPUT WRITERS
# ══════════════════════════════════════════════════════════════

def clean_df(df):
    return df.replace([inf, -inf], np.nan)


def write_haddock_files(haddock_pairs, df_poses, haddock_dir):
    """Writes HADDOCK active/passive residue files for each pair."""
    haddock_dir.mkdir(parents=True, exist_ok=True)

    def passives(resids):
        nbrs = {r + d for r in resids for d in (-1, 1) if r + d > 0}
        return sorted(nbrs - set(resids))

    for pair in haddock_pairs:
        apt, var, bf = pair['aptamer'], pair['variant'], pair['BestPoseFile']
        base = f"{var}__{apt}"
        mask = ((df_poses['variant'] == var) &
                (df_poses['aptamer'] == apt) &
                (df_poses['pose_file'] == bf))
        contacts = []
        row = df_poses[mask]
        if not row.empty:
            try:
                contacts = json.loads(row.iloc[0]['contacts'])
            except Exception:
                pass

        rec_r = sorted({c['phos_resid'] for c in contacts})
        apt_r = sorted({c['apt_resid']  for c in contacts
                        if 'apt_resid' in c})

        for fname, resids in [
            (f"{base}__rec_actives.txt",  rec_r),
            (f"{base}__rec_passives.txt", passives(rec_r)),
            (f"{base}__apt_actives.txt",  apt_r),
            (f"{base}__apt_passives.txt", passives(apt_r)),
        ]:
            with open(haddock_dir / fname, 'w') as f:
                f.writelines(f"{r}\n" for r in resids)


def write_all_outputs(df_poses, df_jobs, df_apt,
                      tier1, tier2, combined, haddock_pairs,
                      hit_matrix, output_prefix):
    ts      = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(f"{output_prefix}_outputs_{ts}")
    out_dir.mkdir(exist_ok=True)

    clean_df(df_poses).to_csv(
        out_dir / f"{output_prefix}.all_poses.csv", index=False)
    clean_df(df_jobs).to_csv(
        out_dir / f"{output_prefix}.per_job_scores.csv", index=False)
    hit_matrix.to_csv(out_dir / "hit_matrix.csv")

    xl = out_dir / f"{output_prefix}.shortlist.xlsx"
    with pd.ExcelWriter(xl, engine='openpyxl') as writer:
        clean_df(df_apt).sort_values(
            'BestHdockScore', ascending=True
        ).to_excel(writer, sheet_name='GlobalRanking', index=False)

        clean_df(tier1).to_excel(
            writer, sheet_name='Tier1_Broad', index=False)
        clean_df(tier2).to_excel(
            writer, sheet_name='Tier2_Specific', index=False)
        clean_df(combined).to_excel(
            writer, sheet_name='Selected_Aptamers', index=False)
        pd.DataFrame(haddock_pairs).to_excel(
            writer, sheet_name='HADDOCK_Pairs', index=False)
        clean_df(df_jobs).sort_values(
            ['variant', 'CSS'], ascending=[True, False]
        ).to_excel(writer, sheet_name='PerJob_Scores', index=False)
        hit_matrix.to_excel(writer, sheet_name='HitMatrix')

    write_haddock_files(
        haddock_pairs, df_poses, out_dir / "haddock_inputs")

    logging.info(f"Outputs written to: {out_dir}")
    return out_dir


# ══════════════════════════════════════════════════════════════
#   FOLDER DISCOVERY
# ══════════════════════════════════════════════════════════════

def discover_jobs(input_dir):
    """Discovers all (variant, aptamer_dir, out_file) triples."""
    jobs = []
    for vdir in sorted(Path(input_dir).iterdir()):
        if not vdir.is_dir():
            continue
        for adir in sorted(vdir.iterdir()):
            if not adir.is_dir():
                continue
            out = vdir / f"{adir.name}_dock.out"
            jobs.append((vdir.name, adir, out if out.exists() else None))
    return jobs


# ══════════════════════════════════════════════════════════════
#   MAIN
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="HDOCK shortlisting — two-tier aptamer-centric ranking")
    parser.add_argument('--input_dir',     required=True)
    parser.add_argument('--phospho_json',  required=True)
    parser.add_argument('--output_prefix', default='wave1_shortlist')
    parser.add_argument('--top_n_poses',   type=int, default=10)
    parser.add_argument('--processes',     type=int,
                        default=max(1, cpu_count() - 2))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s')

    logging.info("=" * 55)
    logging.info("  HDOCK Shortlisting — Two-Tier Aptamer-Centric")
    logging.info("=" * 55)
    logging.info(f"Proximity cutoff : {PROX_CUTOFF} A "
                 f"(Dominguez et al. 2003 JACS)")
    logging.info(f"Direct cutoff    : {DIRECT_CUTOFF} A "
                 f"(Huang & Zou 2014 NAR)")
    logging.info(f"Tier 1           : >= {MIN_HITS_TIER1} variants hit")
    logging.info(f"Tier 2           : top {TOP_PERCENT_TIER2}% HDOCK score")
    logging.info(f"HADDOCK cap      : {MAX_HADDOCK_INPUT} aptamers")

    with open(args.phospho_json, 'r') as f:
        phospho_map = {k: [int(x) for x in v]
                       for k, v in json.load(f).items()}

    jobs = discover_jobs(args.input_dir)
    if not jobs:
        logging.error("No folders found. Check input_dir structure.")
        return
    logging.info(f"Found {len(jobs)} aptamer folders to analyze")

    worker_args = [
        (variant, apt_dir, out_file,
         phospho_map.get(variant, []), args.top_n_poses)
        for variant, apt_dir, out_file in jobs
    ]

    # Parallel pose analysis
    logging.info("Analyzing poses (this is the heavy step)...")
    all_records = []
    if args.processes <= 1:
        for wa in worker_args:
            all_records.extend(analyze_aptamer(wa))
    else:
        with Pool(processes=args.processes) as pool:
            for batch in pool.imap_unordered(analyze_aptamer, worker_args):
                all_records.extend(batch)

    df_poses = pd.DataFrame([r for r in all_records
                             if r.get('error') is None])
    if df_poses.empty:
        logging.error("No valid pose records produced.")
        return
    logging.info(f"Valid pose records: {len(df_poses)}")

    logging.info("Aggregating per job...")
    df_jobs = aggregate_per_job(df_poses)
    df_jobs = compute_css(df_jobs)

    logging.info("Aggregating per aptamer (cross-variant)...")
    df_apt = aggregate_per_aptamer(df_jobs)

    logging.info("Applying two-tier selection...")
    tier1, tier2, combined, haddock_pairs = select_tiers(df_apt, df_jobs)

    hit_matrix = build_hit_matrix(df_jobs)

    out_dir = write_all_outputs(
        df_poses, df_jobs, df_apt,
        tier1, tier2, combined, haddock_pairs,
        hit_matrix, args.output_prefix)

    # ── Summary printout ──────────────────────────────────
    print("\n" + "=" * 60)
    print("  SHORTLISTING COMPLETE")
    print("=" * 60)
    print(f"  Aptamers analyzed        : {len(df_apt)}")
    print(f"  Aptamers with >= 1 hit   : "
          f"{len(df_apt[df_apt['NumHits'] >= 1])}")
    print(f"\n  Tier 1 (>= {MIN_HITS_TIER1} variants)   : {len(tier1)}")
    print(f"  Tier 2 (top {TOP_PERCENT_TIER2}% score)  : {len(tier2)}")
    print(f"  Total selected           : {len(combined)}")
    print(f"  HADDOCK pairs generated  : {len(haddock_pairs)}")
    total_hits   = int((hit_matrix >= 1).values.sum())
    direct_hits  = int((hit_matrix == 2).values.sum())
    print(f"\n  Hit matrix summary:")
    print(f"    Proximal hits (<=6.0A) : {total_hits}")
    print(f"    Direct hits  (<=4.5A)  : {direct_hits}")
    print(f"\n  Hard rejection breakdown:")
    rej = df_jobs[df_jobs['HardReject']]['RejectReason'].value_counts()
    for reason, count in rej.items():
        print(f"    {reason:<30} : {count}")
    if not combined.empty:
        print(f"\n  Top 10 selected aptamers:")
        print(f"  {'Aptamer':<15} {'Class':<6} {'Hits':<6} "
              f"{'MeanHDOCK':<12} {'Tier'}")
        print(f"  {'-'*56}")
        for _, row in combined.head(10).iterrows():
            mh = f"{row['MeanHdockHit']:.2f}" if row['MeanHdockHit'] else 'N/A'
            print(f"  {row['aptamer']:<15} {row['Class']:<6} "
                  f"{row['NumHits']:<6} {mh:<12} {row['Tier']}")
    print(f"\n  Outputs: {out_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
