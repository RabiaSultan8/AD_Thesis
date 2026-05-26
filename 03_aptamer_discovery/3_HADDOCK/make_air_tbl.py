#!/usr/bin/env python3
"""
make_air_tbl.py
Converts 4-file AIR residue lists (.txt) into HADDOCK-format .tbl files.

For each pair in the CSV:
  - rec_actives  → segid A, active
  - rec_passives → segid A, passive
  - apt_actives  → segid B, active
  - apt_passives → segid B, passive

AIR logic (standard HADDOCK):
  Each active residue on A → all (active+passive) on B
  Each active residue on B → all (active+passive) on A

Distance: 2.0 2.0 0.0 (ambig, lower=0, upper=4.0 Å)

Usage:
  python make_air_tbl.py \
      --csv  haddock_pairs_manifest.csv \
      --base /mnt/d/HDOCK \
      --outdir HADDOCK_AIR_tbl
"""

import csv
import argparse
from pathlib import Path


def read_residues(filepath):
    """Read one residue number per line, return sorted list of ints."""
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"AIR file not found: {filepath}")
    residues = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if line:
            residues.append(int(line))
    return sorted(residues)


def make_assign_block(active_res, segid_active,
                      partner_res, segid_partner):
    """
    Build one assign block:
      active residue (segid_active) <-> all partner residues (segid_partner)
    Returns empty string if either list is empty.
    """
    if not active_res or not partner_res:
        return ""

    lines = []
    for act in active_res:
        lines.append(f"assign ( resid {act:>4}  and segid {segid_active} )")
        lines.append( "       (")
        for i, par in enumerate(partner_res):
            prefix = "        " if i == 0 else "     or\n        "
            if i == 0:
                lines.append(f"        ( resid {par:>4}  and segid {segid_partner} )")
            else:
                lines.append(f"     or")
                lines.append(f"        ( resid {par:>4}  and segid {segid_partner} )")
        lines.append( "       )  2.0 2.0 0.0")
        lines.append("!")
    return "\n".join(lines) + "\n"


def make_tbl(rec_act, rec_pas, apt_act, apt_pas):
    """Build the full .tbl content for one pair."""
    rec_all = sorted(set(rec_act + rec_pas))   # active+passive on receptor
    apt_all = sorted(set(apt_act + apt_pas))   # active+passive on aptamer

    tbl = []
    tbl.append("! HADDOCK AIR restraints — receptor actives vs aptamer (active+passive)")
    tbl.append("!")
    block1 = make_assign_block(rec_act, "A", apt_all, "B")
    if block1:
        tbl.append(block1)

    tbl.append("! HADDOCK AIR restraints — aptamer actives vs receptor (active+passive)")
    tbl.append("!")
    block2 = make_assign_block(apt_act, "B", rec_all, "A")
    if block2:
        tbl.append(block2)

    return "\n".join(tbl) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv",    required=True,
                        help="haddock_pairs_manifest.csv")
    parser.add_argument("--base",   default=".",
                        help="Base directory where wave1_shortlist_outputs... lives")
    parser.add_argument("--outdir", default="HADDOCK_AIR_tbl",
                        help="Output directory for .tbl files")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base   = Path(args.base)

    with open(args.csv) as f:
        reader = csv.DictReader(f)
        rows   = list(reader)

    print(f"\n{'Pair':<40} {'RecAct':>7} {'RecPas':>7} {'AptAct':>7} {'AptPas':>7}  Status")
    print("-" * 90)

    ok = 0
    errors = []

    for row in rows:
        aptamer = row["aptamer"]
        variant = row["BestVariant"]
        pair_id = f"{variant}__{aptamer}"

        try:
            rec_act = read_residues(base / row["AIR_rec_act"])
            rec_pas = read_residues(base / row["AIR_rec_pas"])
            apt_act = read_residues(base / row["AIR_apt_act"])
            apt_pas = read_residues(base / row["AIR_apt_pas"])

            tbl_content = make_tbl(rec_act, rec_pas, apt_act, apt_pas)

            out_file = outdir / f"{pair_id}__ambig.tbl"
            out_file.write_text(tbl_content)

            print(f"{pair_id:<40} {len(rec_act):>7} {len(rec_pas):>7} "
                  f"{len(apt_act):>7} {len(apt_pas):>7}  OK -> {out_file.name}")
            ok += 1

        except Exception as e:
            msg = f"ERROR: {e}"
            print(f"{pair_id:<40} {'':>7} {'':>7} {'':>7} {'':>7}  {msg}")
            errors.append((pair_id, str(e)))

    print(f"\n{ok}/{len(rows)} .tbl files written to {outdir}/")
    if errors:
        print("\nFailed pairs:")
        for p, e in errors:
            print(f"  {p}: {e}")
    else:
        print("All pairs OK.")


if __name__ == "__main__":
    main()
