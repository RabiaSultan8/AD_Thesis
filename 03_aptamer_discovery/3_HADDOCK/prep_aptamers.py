#!/usr/bin/env python3
"""
prep_aptamers.py
Prepares aptamer PDB files for HADDOCK docking.

Actions:
  - Renames chain A -> B (col 22)
  - Adds END record if missing
  - Validates: OP1/OP2 naming, O2' presence, single-letter residue names

Usage: python prep_aptamers.py --indir HADDOCK_aptamers --outdir HADDOCK_aptamers/prepared
"""

import os
import argparse
from pathlib import Path

VALID_RNA_RESIDUES = {'A', 'U', 'G', 'C'}

def process_aptamer(in_path, out_path):
    issues = []
    warnings = []
    out_lines = []
    chain_fixed = 0
    has_o2prime = False
    has_op1op2  = False
    bad_residues = set()

    with open(in_path) as f:
        lines = f.readlines()

    for line in lines:
        record = line[:6].strip()

        if record in ("ATOM", "HETATM"):
            line = list(line)
            while len(line) < 80:
                line.append(' ')

            # Fix chain A -> B (index 21)
            if line[21] == 'A':
                line[21] = 'B'
                chain_fixed += 1

            atom_name = ''.join(line[12:16]).strip()
            resname   = ''.join(line[17:21]).strip()

            # Validate residue names
            if resname not in VALID_RNA_RESIDUES:
                bad_residues.add(resname)

            # Check OP1/OP2
            if atom_name in ('OP1', 'OP2'):
                has_op1op2 = True

            # Check O2'
            if atom_name == "O2'":
                has_o2prime = True

            line = ''.join(line)
            if not line.endswith('\n'):
                line += '\n'
            out_lines.append(line)

        elif record in ("TER", "REMARK", "END"):
            out_lines.append(line)

        else:
            out_lines.append(line)

    # Add END if missing
    if not any(l.startswith('END') for l in out_lines):
        out_lines.append('END\n')
        issues.append("Added END")

    # Warnings
    if not has_op1op2:
        warnings.append("WARN: No OP1/OP2 found — check phosphate naming")
    if not has_o2prime:
        warnings.append("WARN: No O2' found — may not be RNA")
    if bad_residues:
        warnings.append(f"WARN: Non-RNA residues found: {bad_residues}")

    with open(out_path, 'w') as f:
        f.writelines(out_lines)

    return {
        "chain_fixed": chain_fixed,
        "has_op1op2":  has_op1op2,
        "has_o2prime": has_o2prime,
        "issues":      issues,
        "warnings":    warnings
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir",  default="HADDOCK_aptamers")
    parser.add_argument("--outdir", default="HADDOCK_aptamers/prepared")
    args = parser.parse_args()

    indir  = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(indir.glob("*.pdb"))
    if not pdbs:
        print(f"No PDB files found in {indir}")
        return

    print(f"\n{'File':<35} {'ChainFix':>9} {'OP1/OP2':>8} {'O2prime':>8}  Notes")
    print("-" * 80)

    all_ok = True
    for pdb in pdbs:
        out_path = outdir / pdb.name
        stats = process_aptamer(pdb, out_path)
        notes = "; ".join(stats["issues"] + stats["warnings"]) or "OK"
        op_flag = "✓" if stats["has_op1op2"]  else "✗"
        o2_flag = "✓" if stats["has_o2prime"] else "✗"
        if stats["warnings"]:
            all_ok = False
        print(f"{pdb.name:<35} {stats['chain_fixed']:>9} {op_flag:>8} {o2_flag:>8}  {notes}")

    status = "All files OK." if all_ok else "Some files have warnings — review before running HADDOCK."
    print(f"\nDone. {len(pdbs)} files written to {outdir}/\n{status}")

if __name__ == "__main__":
    main()
