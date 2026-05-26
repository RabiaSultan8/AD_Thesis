#!/usr/bin/env python3
"""
setup_haddock_runs.py
Creates HADDOCK Phase 1 run directories for all 40 aptamer-pTau pairs.

Directory structure created:
  haddock_runs/
    pTau_2.2__Aptamer_3448/
      pTau_2.2.pdb        (receptor, chain A)
      Aptamer_3448.pdb    (aptamer,  chain B)
      ambig.tbl           (AIR restraints)
      run.param
    pTau_3.3__Aptamer_4700/
      ...

Usage:
  python setup_haddock_runs.py \
      --csv        haddock_pairs_manifest.csv \
      --ptau_dir   HADDOCK_pTau_variants \
      --apt_dir    HADDOCK_aptamers \
      --tbl_dir    HADDOCK_AIR_tbl \
      --haddock    /mnt/d/Installations/haddock2.5-2025-08 \
      --outdir     haddock_runs
"""

import csv
import shutil
import argparse
from pathlib import Path


RUNPARAM_TEMPLATE = """\
AMBIG_TBL=./ambig.tbl
HADDOCK_DIR={haddock_dir}/
N_COMP=2
PDB_FILE1=./{ptau_pdb}
PDB_FILE2=./{apt_pdb}
PROJECT_DIR=./
PROT_SEGID_1=A
PROT_SEGID_2=B
RUN_NUMBER=1
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv",      required=True)
    parser.add_argument("--ptau_dir", default="HADDOCK_pTau_variants")
    parser.add_argument("--apt_dir",  default="HADDOCK_aptamers")
    parser.add_argument("--tbl_dir",  default="HADDOCK_AIR_tbl")
    parser.add_argument("--haddock",  required=True,
                        help="Absolute path to HADDOCK installation")
    parser.add_argument("--outdir",   default="haddock_runs")
    args = parser.parse_args()

    ptau_dir  = Path(args.ptau_dir)
    apt_dir   = Path(args.apt_dir)
    tbl_dir   = Path(args.tbl_dir)
    haddock   = str(Path(args.haddock).resolve())
    outdir    = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.csv) as f:
        rows = list(csv.DictReader(f))

    print(f"\n{'Run Directory':<45} {'Status'}")
    print("-" * 65)

    ok = 0
    errors = []

    for row in rows:
        variant = row["BestVariant"]
        aptamer = row["aptamer"]
        pair_id = f"{variant}__{aptamer}"

        ptau_src = ptau_dir / f"{variant}.pdb"
        apt_src  = apt_dir  / f"{aptamer}.pdb"
        tbl_src  = tbl_dir  / f"{pair_id}__ambig.tbl"

        missing = []
        if not ptau_src.exists(): missing.append(str(ptau_src))
        if not apt_src.exists():  missing.append(str(apt_src))
        if not tbl_src.exists():  missing.append(str(tbl_src))

        if missing:
            msg = "MISSING: " + ", ".join(missing)
            print(f"{pair_id:<45} {msg}")
            errors.append((pair_id, msg))
            continue

        # Create run directory
        run_dir = outdir / pair_id
        run_dir.mkdir(exist_ok=True)

        # Copy files
        shutil.copy2(ptau_src, run_dir / f"{variant}.pdb")
        shutil.copy2(apt_src,  run_dir / f"{aptamer}.pdb")
        shutil.copy2(tbl_src,  run_dir / "ambig.tbl")

        # Write run.param
        run_param = RUNPARAM_TEMPLATE.format(
            haddock_dir = haddock,
            ptau_pdb    = f"{variant}.pdb",
            apt_pdb     = f"{aptamer}.pdb",
        )
        (run_dir / "run.param").write_text(run_param)

        print(f"{pair_id:<45} OK")
        ok += 1

    print(f"\n{ok}/{len(rows)} run directories created in {outdir}/")
    if errors:
        print("\nFailed:")
        for p, e in errors:
            print(f"  {p}: {e}")
    else:
        print("All OK.")

    if ok > 0:
        print(f"""
Next steps:
  1. Patch the master run.cns (once):
     python patch_master_runcns.py --haddock {haddock}

  2. Source HADDOCK environment:
     source {haddock}/haddock_configure.sh

  3. Run a test docking first (single pair):
     cd {outdir}/pTau_2.2__Aptamer_3448
     haddock2.5

  4. If test passes, launch all 40 in parallel (see batch script).
""")


if __name__ == "__main__":
    main()
