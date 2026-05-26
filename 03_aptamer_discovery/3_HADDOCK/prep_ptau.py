#!/usr/bin/env python3
"""
prep_ptau.py
Prepares pTau variant PDB files for HADDOCK docking.

Actions:
  - Strips MODEL, ENDMDL, CRYST1, TITLE, REMARK records
  - Sets chain ID to 'A' (col 22)
  - Renames residue TPO -> TOP
  - Deletes H3T atoms (undefined in HADDOCK toppar)
  - Adds END record if missing

Usage: python prep_ptau.py --indir HADDOCK_pTau_variants --outdir HADDOCK_pTau_variants/prepared
"""

import os
import argparse
from pathlib import Path

STRIP_RECORDS = {"MODEL", "ENDMDL", "CRYST1", "TITLE", "REMARK", "MASTER", "AUTHOR",
                 "EXPDTA", "JRNL", "SEQRES", "SOURCE", "KEYWDS", "REVDAT", "HEADER"}

def process_ptau(in_path, out_path):
    issues = []
    out_lines = []
    tpo_count = 0
    h3t_count = 0
    chain_fixed = 0

    with open(in_path) as f:
        lines = f.readlines()

    for line in lines:
        record = line[:6].strip()

        # Strip unwanted records
        if record in STRIP_RECORDS:
            continue

        # Process ATOM and HETATM
        if record in ("ATOM", "HETATM"):
            line = list(line)

            # Pad line to at least 80 chars
            while len(line) < 80:
                line.append(' ')

            # Fix chain ID (col 22, index 21)
            if line[21] != 'A':
                chain_fixed += 1
                line[21] = 'A'

            # Rename TPO -> TOP (residue name at cols 18-20, indices 17-19)
            resname = ''.join(line[17:21]).strip()
            if resname == 'TPO':
                line[17] = 'T'
                line[18] = 'O'
                line[19] = 'P'
                tpo_count += 1

            # Skip H3T atoms
            atom_name = ''.join(line[12:16]).strip()
            if atom_name == 'H3T':
                h3t_count += 1
                continue

            line = ''.join(line)
            # Ensure newline
            if not line.endswith('\n'):
                line += '\n'
            out_lines.append(line)

        elif record == "TER":
            out_lines.append(line)

        elif record == "END":
            out_lines.append(line)

        else:
            # Keep CONECT, etc. if present
            out_lines.append(line)

    # Add END if missing
    if not any(l.startswith('END') for l in out_lines):
        out_lines.append('END\n')
        issues.append("Added missing END record")

    with open(out_path, 'w') as f:
        f.writelines(out_lines)

    return {
        "tpo_renamed": tpo_count,
        "h3t_removed": h3t_count,
        "chain_fixed": chain_fixed,
        "issues": issues
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir",  default="healthy_tau")
    parser.add_argument("--outdir", default="healthy_tau/prepared")
    args = parser.parse_args()

    indir  = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(indir.glob("*.pdb"))
    if not pdbs:
        print(f"No PDB files found in {indir}")
        return

    print(f"\n{'File':<30} {'Chain':>6} {'TPO→TOP':>8} {'H3T del':>8} {'Notes'}")
    print("-" * 75)

    for pdb in pdbs:
        out_path = outdir / pdb.name
        stats = process_ptau(pdb, out_path)
        notes = "; ".join(stats["issues"]) if stats["issues"] else "OK"
        print(f"{pdb.name:<30} {stats['chain_fixed']:>6} {stats['tpo_renamed']:>8} "
              f"{stats['h3t_removed']:>8}   {notes}")

    print(f"\nDone. {len(pdbs)} files written to {outdir}/")

if __name__ == "__main__":
    main()
