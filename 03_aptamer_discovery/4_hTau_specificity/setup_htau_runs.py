#!/usr/bin/env python3
"""
setup_htau_runs.py
Creates HADDOCK run directories for all 20 hTau control dockings.
"""

import shutil
from pathlib import Path

BASE     = Path("/home/zxkj006/haddock_work/hTau_dockings")
HTAU_DIR = BASE / "healthy_tau"
APT_DIR  = BASE / "top20_aptamers"
TBL_DIR  = BASE / "hTau_dockings_AIR_tbl"
RUNS_DIR = BASE / "htau_runs"
HADDOCK  = "/home/zxkj006/haddock_work/haddock2.5-2025-08"

RUNS_DIR.mkdir(exist_ok=True)

RUNPARAM_TEMPLATE = """\
AMBIG_TBL=./ambig.tbl
HADDOCK_DIR={haddock_dir}/
N_COMP=2
PDB_FILE1=./{htau_pdb}
PDB_FILE2=./{apt_pdb}
PROJECT_DIR=./
PROT_SEGID_1=A
PROT_SEGID_2=B
RUN_NUMBER=1
"""

top20 = [
    ("pTau_2.9", "Aptamer_1152"),
    ("pTau_3.3", "Aptamer_2651"),
    ("pTau_2.9", "Aptamer_184"),
    ("pTau_1.8", "Aptamer_1230"),
    ("pTau_3.8", "Aptamer_3305"),
    ("pTau_4.8", "Aptamer_2414"),
    ("pTau_2.4", "Aptamer_138"),
    ("pTau_2.4", "Aptamer_15"),
    ("pTau_2.8", "Aptamer_1223"),
    ("pTau_4.3", "Aptamer_212"),
    ("pTau_3.3", "Aptamer_105"),
    ("pTau_2.4", "Aptamer_80"),
    ("pTau_4.8", "Aptamer_4154"),
    ("pTau_2.9", "Aptamer_838"),
    ("pTau_2.8", "Aptamer_2711"),
    ("pTau_2.9", "Aptamer_717"),
    ("pTau_2.9", "Aptamer_1743"),
    ("pTau_2.2", "Aptamer_1210"),
    ("pTau_1.10","Aptamer_164"),
    ("pTau_2.1", "Aptamer_3553"),
]

print(f"\n{'Pair':<40} {'Status'}")
print("-" * 55)

ok = 0
for ptau, apt in top20:
    conformer = ptau.split("_")[1].split(".")[0]
    htau_name = f"hTau_{conformer}"
    pair_id   = f"{htau_name}__{apt}"

    htau_src = HTAU_DIR / f"{htau_name}.pdb"
    apt_src  = APT_DIR  / f"{apt}.pdb"
    tbl_src  = TBL_DIR  / f"{pair_id}__ambig.tbl"

    run_dir = RUNS_DIR / pair_id
    run_dir.mkdir(exist_ok=True)

    shutil.copy2(htau_src, run_dir / f"{htau_name}.pdb")
    shutil.copy2(apt_src,  run_dir / f"{apt}.pdb")
    shutil.copy2(tbl_src,  run_dir / "ambig.tbl")

    run_param = RUNPARAM_TEMPLATE.format(
        haddock_dir = HADDOCK,
        htau_pdb    = f"{htau_name}.pdb",
        apt_pdb     = f"{apt}.pdb",
    )
    (run_dir / "run.param").write_text(run_param)

    print(f"{pair_id:<40} OK")
    ok += 1

print(f"\n{ok}/20 run directories created in {RUNS_DIR}/")
