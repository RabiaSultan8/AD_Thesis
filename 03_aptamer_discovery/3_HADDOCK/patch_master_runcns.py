#!/usr/bin/env python3
"""
patch_master_runcns.py
Patches protocols/run.cns in the HADDOCK installation for RNA docking.

Run ONCE after installation. All subsequent runs will inherit these settings.

Changes applied:
  dna_mol2         false  -> true
  prot_top_mol2    protein -> dna-rna-allatom-hj-opls-1.3.top
  prot_link_mol2   protein -> dna-rna-1.3.link
  prot_par_mol2    protein -> dna-rna-allatom-hj-opls-1.3.param
  epsilon_0        10.0   -> 78.0
  epsilon_1        1.0    -> 78.0
  dielec_0         rdie   -> cdie
  dielec_1         rdie   -> cdie
  dnarest_on       false  -> true
  noecv            true   -> false

Usage:
  python patch_master_runcns.py --haddock /mnt/d/Installations/haddock2.5-2025-08
"""

import re
import shutil
import argparse
from pathlib import Path
from datetime import datetime

# Each entry: (description, regex_pattern, replacement_string)
PATCHES = [
    (
        "dna_mol2: false -> true",
        r'(\{===>}\s*dna_mol2=)false(;)',
        r'\1true\2'
    ),
    (
        "prot_top_mol2: protein -> RNA topology",
        r'(\{===>}\s*prot_top_mol2=)"protein-allhdg5-4\.top"(;)',
        r'\1"dna-rna-allatom-hj-opls-1.3.top"\2'
    ),
    (
        "prot_link_mol2: protein -> RNA linkage",
        r'(\{===>}\s*prot_link_mol2=)"protein-allhdg5-4-noter\.link"(;)',
        r'\1"dna-rna-1.3.link"\2'
    ),
    (
        "prot_par_mol2: protein -> RNA parameters",
        r'(\{===>}\s*prot_par_mol2=)"protein-allhdg5-4\.param"(;)',
        r'\1"dna-rna-allatom-hj-opls-1.3.param"\2'
    ),
    (
        "epsilon_0: 10.0 -> 78.0",
        r'(\{===>}\s*epsilon_0=)10\.0(;)',
        r'\g<1>78.0\2'
    ),
    (
        "epsilon_1: 1.0 -> 78.0",
        r'(\{===>}\s*epsilon_1=)1\.0(;)',
        r'\g<1>78.0\2'
    ),
    (
        "dielec_0: rdie -> cdie",
        r'(\{===>}\s*dielec_0=)rdie(;)',
        r'\1cdie\2'
    ),
    (
        "dielec_1: rdie -> cdie",
        r'(\{===>}\s*dielec_1=)rdie(;)',
        r'\1cdie\2'
    ),
    (
        "dnarest_on: false -> true",
        r'(\{===>}\s*dnarest_on=)false(;)',
        r'\1true\2'
    ),
    (
        "noecv: true -> false",
        r'(\{===>}\s*noecv=)true(;)',
        r'\1false\2'
    ),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--haddock", required=True,
                        help="Path to HADDOCK installation directory")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would change without writing")
    args = parser.parse_args()

    runcns = Path(args.haddock) / "protocols" / "run.cns"
    if not runcns.exists():
        print(f"ERROR: run.cns not found at {runcns}")
        return

    # Backup
    backup = runcns.with_suffix(f".cns.bak_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    if not args.dry_run:
        shutil.copy2(runcns, backup)
        print(f"Backup: {backup}")

    content = runcns.read_text()
    original = content

    print(f"\nPatching: {runcns}\n")
    print(f"{'Change':<50} {'Status':>10}")
    print("-" * 62)

    for desc, pattern, replacement in PATCHES:
        new_content, n = re.subn(pattern, replacement, content, flags=re.MULTILINE)
        status = f"OK ({n} match)" if n > 0 else "NOT FOUND"
        print(f"{desc:<50} {status:>10}")
        if n == 0:
            print(f"  WARNING: Pattern not matched — check run.cns manually")
        content = new_content

    if content == original:
        print("\nNo changes made — file may already be patched.")
        return

    if args.dry_run:
        print("\n[DRY RUN] No file written.")
    else:
        runcns.write_text(content)
        print(f"\nrun.cns patched successfully.")


if __name__ == "__main__":
    main()
