# =============================================================
#   RNA APTAMER COMPUTATIONAL DESIGN PIPELINE
#   Target: Phosphorylated Tau Protein (Alzheimer's Disease)
#   Scripts: 1 → 2 → 3 → 4 → 5 → 6 → 7
#   Run in sequence. Each script reads outputs of the previous.
# =============================================================
#   PIPELINE OVERVIEW
#   Script 1  →  Generate 5,000-sequence aptamer library
#   Script 2  →  Run RNAfold on all sequences
#   Script 3  →  Apply six biophysical filters
#   Script 4  →  CD-HIT-EST 80% redundancy removal
#   Script 5  →  RNApdist ensemble clustering
#   Script 6  →  Compile final annotated library
#   Script 7  →  Generate all publication figures
# =============================================================


# ─────────────────────────────────────────────────────────────
#   SCRIPT 1 — RNA Aptamer Library Generator
# ─────────────────────────────────────────────────────────────
"""
Generates a 5,000-sequence RNA aptamer library for
computational screening against phosphorylated Tau protein.

Target phosphorylation sites:
    Early-stage markers  : pT181, pT217, pT231, pS202
    Microtubule-binding  : pS262, pS305, pS356
    Aggregation drivers  : pS396, pS404, pS409, pS412, pS422

Library composition (1,250 sequences per class):
    Class A (Aptamer_1    – Aptamer_1250) : G4-Core loop
    Class B (Aptamer_1251 – Aptamer_2500) : UGGG-Anchor loop
    Class C (Aptamer_2501 – Aptamer_3750) : Distributed-G loop
    Class D (Aptamer_3751 – Aptamer_5000) : Mixed scaffold loop

Output:
    aptamer_library.fasta
"""

import random
from collections import defaultdict

# ── Parameters ───────────────────────────────────────────────

RANDOM_SEED  = 1122
NUM_TOTAL    = 5000
CLASS_SIZE   = 1250
LENGTHS      = [28, 32, 36, 40]

# Stem length per total length; loop = total − (2 × stem) ∈ [10, 14]
STEM_LENGTH_MAP = {28: 8, 32: 9, 36: 11, 40: 13}

STEM_GC_MIN  = 0.60
STEM_GC_MAX  = 0.70
STEM_GC_POOL = ['G', 'C']
STEM_AU_POOL = ['A', 'U']

LOOP_G_MIN   = 0.40
LOOP_G_MAX   = 0.55
LOOP_NON_G   = ['A', 'U', 'C']

MAX_POLY_G   = 5
MAX_POLY_C   = 4

OUTPUT_FASTA = "aptamer_library.fasta"

random.seed(RANDOM_SEED)

# ── Helper functions ─────────────────────────────────────────

def wc_complement(base):
    """Returns the Watson-Crick RNA complement of a base."""
    return {'G': 'C', 'C': 'G', 'A': 'U', 'U': 'A'}.get(base)


def reverse_complement(seq):
    """Returns the RNA reverse complement of a stem sequence."""
    rc = [wc_complement(b) for b in reversed(seq)]
    return None if None in rc else ''.join(rc)


def generate_stem(length):
    """Generates a random 5' stem with GC content in [60%, 70%]."""
    gc_ratio = random.uniform(STEM_GC_MIN, STEM_GC_MAX)
    num_gc   = round(length * gc_ratio)
    bases    = (random.choices(STEM_GC_POOL, k=num_gc) +
                random.choices(STEM_AU_POOL, k=length - num_gc))
    random.shuffle(bases)
    return ''.join(bases)


def apply_gu_wobble(three_stem):
    """
    Introduces one G-U wobble by replacing one internal C → U.
    Terminal positions are left unchanged. Used only for Class C.
    """
    stem = list(three_stem)
    eligible = [i for i in range(1, len(stem) - 1) if stem[i] == 'C']
    if eligible:
        stem[random.choice(eligible)] = 'U'
    return ''.join(stem)


def check_poly_runs(seq):
    """Returns True if sequence has no poly-G > 5 or poly-C > 4."""
    return ('G' * (MAX_POLY_G + 1) not in seq and
            'C' * (MAX_POLY_C + 1) not in seq)


def check_loop_g_content(loop):
    """Returns True if loop G-content is within [40%, 55%]."""
    if not loop:
        return False
    return LOOP_G_MIN <= loop.count('G') / len(loop) <= LOOP_G_MAX


def g_enriched_fill(n):
    """Generates n random bases with ~45% G content for loop padding."""
    bases = []
    for _ in range(n):
        r = random.random()
        if   r < 0.45: bases.append('G')
        elif r < 0.65: bases.append('A')
        elif r < 0.82: bases.append('U')
        else:          bases.append('C')
    return ''.join(bases)


# ── Class-specific generators ────────────────────────────────

def make_class_A(length):
    """
    Class A — G4-Core Loop
    Loop contains 2–3 GGG tracts separated by 1–3 non-G spacers,
    forming a G-quadruplex-like motif that coordinates phosphate
    oxygens of phospho-Tau residues.
    """
    stem_len = STEM_LENGTH_MAP[length]
    loop_len = length - 2 * stem_len
    if not (10 <= loop_len <= 14):
        return None

    num_tracts = random.choice([2, 3])
    loop = []
    for i in range(num_tracts):
        loop.extend(['G', 'G', 'G'])
        if i < num_tracts - 1:
            loop.extend(random.choices(LOOP_NON_G, k=random.randint(1, 3)))

    diff = loop_len - len(loop)
    if diff < 0:
        non_g = [i for i, b in enumerate(loop) if b != 'G']
        while len(loop) > loop_len and non_g:
            loop.pop(non_g.pop())
            non_g = [i for i, b in enumerate(loop) if b != 'G']
        if len(loop) != loop_len:
            return None
    elif diff > 0:
        pos = random.randint(0, len(loop))
        loop = loop[:pos] + list(g_enriched_fill(diff)) + loop[pos:]

    loop_seq = ''.join(loop)
    if len(loop_seq) != loop_len or not check_loop_g_content(loop_seq):
        return None

    five_stem  = generate_stem(stem_len)
    three_stem = reverse_complement(five_stem)
    if three_stem is None:
        return None

    seq = five_stem + loop_seq + three_stem
    if not check_poly_runs(seq):
        return None
    return seq, five_stem, three_stem, loop_seq


def make_class_B(length):
    """
    Class B — UGGG-Anchor Loop
    Loop contains a UGGG or UGGAGG motif at a random internal
    position; these motifs make direct H-bond contacts with
    phosphate oxygens through their guanosine residues.
    """
    stem_len  = STEM_LENGTH_MAP[length]
    loop_len  = length - 2 * stem_len
    if not (10 <= loop_len <= 14):
        return None

    motif    = random.choice(['UGGG', 'UGGAGG'])
    remaining = loop_len - len(motif)
    if remaining < 2:
        return None

    prefix_len = random.randint(1, remaining - 1)
    loop_seq   = (g_enriched_fill(prefix_len) + motif +
                  g_enriched_fill(remaining - prefix_len))

    if len(loop_seq) != loop_len or not check_loop_g_content(loop_seq):
        return None

    five_stem  = generate_stem(stem_len)
    three_stem = reverse_complement(five_stem)
    if three_stem is None:
        return None

    seq = five_stem + loop_seq + three_stem
    if not check_poly_runs(seq):
        return None
    return seq, five_stem, three_stem, loop_seq


def make_class_C(length):
    """
    Class C — Distributed-G Loop with G-U Wobble
    G bases are spread evenly across the loop (not clustered),
    creating a broad electrostatic surface. The 3' stem carries
    exactly one internal G-U wobble that promotes tertiary contact
    with disordered Tau loop regions.
    """
    stem_len = STEM_LENGTH_MAP[length]
    loop_len = length - 2 * stem_len
    if not (10 <= loop_len <= 14):
        return None

    g_count   = round(loop_len * random.uniform(LOOP_G_MIN, LOOP_G_MAX))
    g_pos_set = set(random.sample(range(loop_len), g_count))
    loop_seq  = ''.join(
        'G' if i in g_pos_set else random.choice(LOOP_NON_G)
        for i in range(loop_len)
    )
    if not check_loop_g_content(loop_seq):
        return None

    five_stem  = generate_stem(stem_len)
    three_stem = apply_gu_wobble(reverse_complement(five_stem) or '')
    if not three_stem:
        return None

    seq = five_stem + loop_seq + three_stem
    if not check_poly_runs(seq):
        return None
    return seq, five_stem, three_stem, loop_seq


def make_class_D(length):
    """
    Class D — Mixed Scaffold Loop
    Combines GGG (G4 sub-motif) and UGGG (anchor motif) in one
    loop separated by a 2-base non-G linker. This dual-motif
    design covers both phosphate coordination modes simultaneously.
    Core structure: GGG-[2 non-G]-UGGG (10 bases).
    """
    stem_len = STEM_LENGTH_MAP[length]
    loop_len = length - 2 * stem_len
    if not (10 <= loop_len <= 14):
        return None

    linker   = ''.join(random.choices(LOOP_NON_G, k=2))
    core     = 'GGG' + linker + 'UGGG'
    if len(core) > loop_len:
        return None

    fill     = g_enriched_fill(loop_len - len(core))
    loop_seq = fill + core if random.random() < 0.5 else core + fill

    if len(loop_seq) != loop_len or not check_loop_g_content(loop_seq):
        return None

    five_stem  = generate_stem(stem_len)
    three_stem = reverse_complement(five_stem)
    if three_stem is None:
        return None

    seq = five_stem + loop_seq + three_stem
    if not check_poly_runs(seq):
        return None
    return seq, five_stem, three_stem, loop_seq


# ── Library generation ───────────────────────────────────────

GENERATORS   = {'A': make_class_A, 'B': make_class_B,
                'C': make_class_C, 'D': make_class_D}
CLASS_ORDER  = ['A', 'B', 'C', 'D']
CLASS_RANGES = {'A': (1, 1250), 'B': (1251, 2500),
                'C': (2501, 3750), 'D': (3751, 5000)}


def generate_library():
    """
    Generates the full 5,000-sequence library class by class
    (A → B → C → D), each contributing exactly 1,250 unique
    sequences numbered globally Aptamer_1 to Aptamer_5000.
    """
    library    = {}
    rejections = defaultdict(int)

    print("=" * 55)
    print("  RNA Aptamer Library Generator")
    print("  Seed: 1122  |  Target: 5,000 sequences")
    print("  Classes: A (1–1250)  B (1251–2500)")
    print("           C (2501–3750)  D (3751–5000)")
    print("=" * 55)

    for cls in CLASS_ORDER:
        start, end   = CLASS_RANGES[cls]
        target       = end - start + 1
        class_seqs   = set()
        counter      = start
        attempts     = 0

        print(f"\nGenerating Class {cls} (Aptamer_{start} to Aptamer_{end})...")

        while len(class_seqs) < target and attempts < target * 50:
            attempts += 1
            result    = GENERATORS[cls](random.choice(LENGTHS))

            if result is None:
                rejections[f'{cls}_invalid'] += 1
                continue

            seq = result[0]
            if seq in class_seqs or seq in library.values():
                rejections[f'{cls}_duplicate'] += 1
                continue

            apt_id          = f"Aptamer_{counter}"
            library[apt_id] = seq
            class_seqs.add(seq)
            counter += 1

        print(f"  Generated  : {len(class_seqs)}")
        print(f"  Attempts   : {attempts}")
        print(f"  Rejections : "
              f"{rejections.get(f'{cls}_invalid', 0)} invalid, "
              f"{rejections.get(f'{cls}_duplicate', 0)} duplicate")

    print(f"\n{'=' * 55}")
    print(f"  Total sequences generated : {len(library)}")
    print(f"{'=' * 55}")
    return library


def write_fasta(library, filename=OUTPUT_FASTA):
    """Writes library to FASTA file in numerical order."""
    sorted_ids = sorted(library.keys(),
                        key=lambda x: int(x.split('_')[1]))
    with open(filename, 'w') as f:
        for apt_id in sorted_ids:
            f.write(f">{apt_id}\n{library[apt_id]}\n")
    print(f"\nLibrary written to : {filename}")
    print(f"Total sequences    : {len(library)}")


if __name__ == "__main__":
    library = generate_library()
    write_fasta(library)




# ─────────────────────────────────────────────────────────────
#   SCRIPT 2 — RNAfold Batch Runner
# ─────────────────────────────────────────────────────────────
"""
Runs RNAfold on all 5,000 sequences in aptamer_library.fasta
and collects secondary structure predictions.

For each sequence, RNAfold produces:
    - MFE structure     : dot-bracket notation
    - MFE energy        : minimum free energy (kcal/mol)
    - Ensemble energy   : partition function energy
    - Ensemble diversity: structural heterogeneity score

Requirement:
    ViennaRNA installed — verify with: RNAfold --version

Input  : aptamer_library.fasta
Output : rnafold_output.txt
"""

import subprocess
from pathlib import Path

# ── Parameters ───────────────────────────────────────────────

INPUT_FASTA  = "aptamer_library.fasta"
OUTPUT_FILE  = "rnafold_output.txt"

# -p    : compute partition function and ensemble diversity
# -d2   : dangling ends on both sides (standard for short RNAs)
# --noPS: skip PostScript structure diagrams
# --noLP: disallow lonely base pairs
RNAFOLD_CMD  = ["RNAfold", "-p", "-d2", "--noPS", "--noLP"]

# ── Helper functions ─────────────────────────────────────────

def read_fasta(filename):
    """Reads FASTA file; returns ordered list of (header, seq) tuples."""
    records, header, seq = [], None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.upper())
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def run_rnafold(header, sequence):
    """
    Runs RNAfold on a single sequence via stdin.
    Returns raw stdout string, or None on failure.
    """
    try:
        result = subprocess.run(
            RNAFOLD_CMD,
            input=f">{header}\n{sequence}\n",
            capture_output=True, text=True, timeout=30
        )
        return result.stdout if result.returncode == 0 else None
    except subprocess.TimeoutExpired:
        return None
    except FileNotFoundError:
        print("ERROR: RNAfold not found. Install ViennaRNA.")
        raise


# ── Main runner ──────────────────────────────────────────────

def run_batch():
    """
    Processes all sequences sequentially.
    Writes all RNAfold output blocks to OUTPUT_FILE.
    Prints progress every 500 sequences.
    """
    records = read_fasta(INPUT_FASTA)
    total   = len(records)
    failed  = []

    print("=" * 55)
    print("  RNAfold Batch Runner")
    print(f"  Input  : {INPUT_FASTA} ({total} sequences)")
    print(f"  Output : {OUTPUT_FILE}")
    print("=" * 55)

    with open(OUTPUT_FILE, 'w') as out:
        for i, (header, seq) in enumerate(records, 1):
            if i % 500 == 0 or i == 1 or i == total:
                print(f"  Processing {i}/{total} — {header}")

            output = run_rnafold(header, seq)

            if output is None:
                print(f"  WARNING: RNAfold failed for {header} — skipping")
                failed.append(header)
                continue

            out.write(output + "\n")

    print(f"\n{'=' * 55}")
    print(f"  Completed  : {total - len(failed)}/{total}")
    print(f"  Failed     : {len(failed)}")
    if failed:
        print(f"  Failed IDs : {', '.join(failed[:10])}"
              f"{'...' if len(failed) > 10 else ''}")
    print(f"  Output     : {OUTPUT_FILE}")
    print(f"{'=' * 55}")


if __name__ == "__main__":
    run_batch()




# ─────────────────────────────────────────────────────────────
#   SCRIPT 3 — Multi-Criterion Biophysical Filter
# ─────────────────────────────────────────────────────────────
"""
Parses rnafold_output.txt and applies six sequential filters
to select structurally sound aptamer candidates.

Filter criteria (all six must pass):
    F1 — MFE total         : <= -8.0 kcal/mol
    F2 — MFE per nt        : <= -0.30 kcal/mol/nt
    F3 — Ensemble diversity: <= 8.0
    F4 — Loop size         : 10–14 nt (confirmed from structure)
    F5 — GGG in loop       : >= 1 GGG tract in loop
    F6 — Stem pairs        : >= 6 confirmed base pairs

Inputs  : rnafold_output.txt, aptamer_library.fasta
Outputs : aptamer_filtered.fasta, filter_report.tsv
"""

import re
from collections import defaultdict

# ── Parameters ───────────────────────────────────────────────

INPUT_RNAFOLD = "rnafold_output.txt"
INPUT_FASTA   = "aptamer_library.fasta"
OUTPUT_FASTA  = "aptamer_filtered.fasta"
OUTPUT_REPORT = "filter_report.tsv"

F1_MFE_TOTAL      = -8.0
F2_MFE_PER_NT     = -0.30
F3_DIVERSITY_MAX  =  8.0
F4_LOOP_MIN       =  10
F4_LOOP_MAX       =  14
F5_GGG_REQUIRED   =  1
F6_STEM_PAIRS_MIN =  6

# ── Helper functions ─────────────────────────────────────────

def read_fasta(filename):
    """Reads FASTA file into a dict: header → sequence."""
    seqs, header = {}, None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                seqs[header] = ''
            elif header:
                seqs[header] += line.upper()
    return seqs


def parse_dot_bracket(structure):
    """
    Extracts loop and stem info from a dot-bracket string.
    The main loop is the longest dot run enclosed by ( ).
    Returns dict with loop_start, loop_end, loop_len, stem_pairs,
    or None if no valid loop is found.
    """
    stem_pairs = structure.count('(')
    dot_runs   = [(m.start(), m.end() - 1, m.end() - m.start())
                  for m in re.finditer(r'\.+', structure)]
    if not dot_runs:
        return None

    valid_loops = [
        (s, e, l) for s, e, l in dot_runs
        if s > 0 and e < len(structure) - 1
        and structure[s - 1] == '(' and structure[e + 1] == ')'
    ]
    if not valid_loops:
        return None

    loop_start, loop_end, loop_len = max(valid_loops, key=lambda x: x[2])
    return {'loop_start': loop_start, 'loop_end': loop_end,
            'loop_len': loop_len, 'stem_pairs': stem_pairs}


def count_ggg_in_loop(sequence, loop_start, loop_end):
    """Counts GGG tracts within the loop region."""
    return len(re.findall(r'GGG', sequence[loop_start:loop_end + 1]))


# ── RNAfold output parser ────────────────────────────────────

def parse_rnafold_output(filename):
    """
    Parses rnafold_output.txt into a dict: header → record.
    Each record: {mfe_structure, mfe_energy, diversity}.
    """
    records = {}
    with open(filename, 'r') as f:
        lines = [l.rstrip() for l in f]

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header    = lines[i][1:].strip()
            mfe_struct = mfe_energy = diversity = None
            i += 2   # skip sequence line

            while i < len(lines) and not lines[i].startswith('>'):
                line = lines[i].strip()
                if (mfe_struct is None and
                        re.match(r'^[\(\)\.\{\}\,\|]+\s+\([-\d\.]+\)', line)):
                    parts      = line.split()
                    mfe_struct = parts[0]
                    mfe_energy = float(parts[-1].strip('()'))
                elif 'ensemble diversity' in line.lower():
                    m = re.search(r'ensemble diversity\s+([\d\.]+)',
                                  line, re.IGNORECASE)
                    if m:
                        diversity = float(m.group(1))
                i += 1

            if mfe_struct and mfe_energy is not None and diversity is not None:
                records[header] = {'mfe_structure': mfe_struct,
                                   'mfe_energy'   : mfe_energy,
                                   'diversity'    : diversity}
        else:
            i += 1

    return records


# ── Filter logic ─────────────────────────────────────────────

def apply_filters(sequences, rnafold_records):
    """
    Applies all six filters. Returns (passed list, report list,
    fail_counts dict).
    """
    passed, report = [], []
    fail_counts    = defaultdict(int)

    for header, seq in sequences.items():
        if header not in rnafold_records:
            fail_counts['no_rnafold_record'] += 1
            continue

        rec       = rnafold_records[header]
        structure = rec['mfe_structure']
        mfe       = rec['mfe_energy']
        diversity = rec['diversity']
        length    = len(seq)
        parsed    = parse_dot_bracket(structure)

        row = {
            'aptamer_id' : header,
            'length'     : length,
            'mfe_total'  : mfe,
            'mfe_per_nt' : round(mfe / length, 4),
            'diversity'  : diversity,
            'loop_len'   : parsed['loop_len']   if parsed else 'NA',
            'stem_pairs' : parsed['stem_pairs'] if parsed else 'NA',
            'ggg_in_loop': 0,
            'pass_F1': False, 'pass_F2': False, 'pass_F3': False,
            'pass_F4': False, 'pass_F5': False, 'pass_F6': False,
            'pass_all': False,
        }

        if parsed is None:
            fail_counts['no_valid_loop'] += 1
            report.append(row)
            continue

        ggg = count_ggg_in_loop(seq, parsed['loop_start'], parsed['loop_end'])
        row['ggg_in_loop'] = ggg

        f1 = mfe <= F1_MFE_TOTAL
        f2 = (mfe / length) <= F2_MFE_PER_NT
        f3 = diversity <= F3_DIVERSITY_MAX
        f4 = F4_LOOP_MIN <= parsed['loop_len'] <= F4_LOOP_MAX
        f5 = ggg >= F5_GGG_REQUIRED
        f6 = parsed['stem_pairs'] >= F6_STEM_PAIRS_MIN

        row.update({'pass_F1': f1, 'pass_F2': f2, 'pass_F3': f3,
                    'pass_F4': f4, 'pass_F5': f5, 'pass_F6': f6,
                    'pass_all': all([f1, f2, f3, f4, f5, f6])})

        for tag, flag in [('F1_mfe_total', f1), ('F2_mfe_per_nt', f2),
                          ('F3_diversity', f3), ('F4_loop_size', f4),
                          ('F5_ggg_loop', f5), ('F6_stem_pairs', f6)]:
            if not flag:
                fail_counts[tag] += 1

        if all([f1, f2, f3, f4, f5, f6]):
            passed.append(header)

        report.append(row)

    return passed, report, fail_counts


# ── Output writers ───────────────────────────────────────────

def write_filtered_fasta(passed, sequences, filename=OUTPUT_FASTA):
    """Writes passing sequences to FASTA in numerical order."""
    sorted_passed = sorted(passed, key=lambda x: int(x.split('_')[1]))
    with open(filename, 'w') as f:
        for h in sorted_passed:
            f.write(f">{h}\n{sequences[h]}\n")
    print(f"\nFiltered FASTA  : {filename}")
    print(f"Passing sequences: {len(sorted_passed)}")


def write_report(report, filename=OUTPUT_REPORT):
    """Writes per-aptamer filter metrics to a TSV file."""
    cols = ['aptamer_id', 'length', 'mfe_total', 'mfe_per_nt',
            'diversity', 'loop_len', 'stem_pairs', 'ggg_in_loop',
            'pass_F1', 'pass_F2', 'pass_F3', 'pass_F4',
            'pass_F5', 'pass_F6', 'pass_all']
    with open(filename, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in report:
            f.write('\t'.join(str(row.get(c, 'NA')) for c in cols) + '\n')
    print(f"Filter report   : {filename}")


if __name__ == "__main__":
    print("=" * 55)
    print("  Biophysical Filter — Script 3")
    print("=" * 55)

    sequences       = read_fasta(INPUT_FASTA)
    rnafold_records = parse_rnafold_output(INPUT_RNAFOLD)
    print(f"\n  Sequences loaded  : {len(sequences)}")
    print(f"  RNAfold records   : {len(rnafold_records)}")

    passed, report, fail_counts = apply_filters(sequences, rnafold_records)

    print(f"\n{'=' * 55}")
    print(f"  Input sequences       : {len(sequences)}")
    print(f"  Passed all filters    : {len(passed)}")
    print(f"  Failed                : {len(sequences) - len(passed)}")
    print(f"\n  Failure breakdown:")
    for reason, count in sorted(fail_counts.items(), key=lambda x: -x[1]):
        print(f"    {reason:<25} : {count}")
    print(f"{'=' * 55}")

    write_filtered_fasta(passed, sequences)
    write_report(report)




# ─────────────────────────────────────────────────────────────
#   SCRIPT 4 — CD-HIT-EST Sequence Redundancy Removal
# ─────────────────────────────────────────────────────────────
"""
Runs CD-HIT-EST on aptamer_filtered.fasta to remove sequences
sharing >= 80% identity, retaining one representative per cluster.

This eliminates near-identical sequences before the more expensive
RNApdist ensemble clustering in Script 5.

Requirement:
    CD-HIT installed — verify with: cd-hit-est -h
    Install with   : conda install -c bioconda cd-hit

Input  : aptamer_filtered.fasta   (1,179 sequences from Script 3)
Outputs: aptamer_cdhit.fasta
         aptamer_cdhit.fasta.clstr
         cdhit_report.txt
"""

import subprocess
import re

# ── Parameters ───────────────────────────────────────────────

INPUT_FASTA        = "aptamer_filtered.fasta"
OUTPUT_PREFIX      = "aptamer_cdhit"
OUTPUT_FASTA       = f"{OUTPUT_PREFIX}.fasta"
CLUSTER_FILE       = f"{OUTPUT_PREFIX}.fasta.clstr"
REPORT_FILE        = "cdhit_report.txt"

IDENTITY_THRESHOLD = 0.80
WORD_SIZE          = 5
THREADS            = 4
MEMORY             = 2000

# ── Helper functions ─────────────────────────────────────────

def read_fasta(filename):
    """Reads FASTA into ordered list of (header, seq) tuples."""
    records, header, seq = [], None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.upper())
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def parse_clstr_file(clstr_file):
    """
    Parses CD-HIT .clstr file.
    Returns (clusters dict, representatives dict, total cluster count).
    """
    clusters, representatives, current = {}, {}, None
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current = int(line.split()[1])
                clusters[current] = []
            elif current is not None and line:
                m = re.search(r'>(\S+?)\.\.\.', line)
                if m:
                    seq_id = m.group(1)
                    clusters[current].append(seq_id)
                    if line.endswith('*'):
                        representatives[current] = seq_id
    return clusters, representatives, len(clusters)


def run_cdhit(input_fasta, output_prefix):
    """Runs CD-HIT-EST. Returns True on success, False on failure."""
    cmd = [
        'cd-hit-est',
        '-i', input_fasta,
        '-o', output_prefix + '.fasta',
        '-c', str(IDENTITY_THRESHOLD),
        '-n', str(WORD_SIZE),
        '-T', str(THREADS),
        '-M', str(MEMORY),
        '-d', '0',
        '-r', '0',
    ]
    print(f"  Running: {' '.join(cmd)}\n")
    try:
        result = subprocess.run(cmd, capture_output=True,
                                text=True, timeout=300)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr}")
            return False
        return True
    except FileNotFoundError:
        print("  ERROR: cd-hit-est not found.")
        return False
    except subprocess.TimeoutExpired:
        print("  ERROR: CD-HIT-EST timed out.")
        return False


def write_report(clusters, input_count, output_count,
                 filename=REPORT_FILE):
    """Writes a human-readable CD-HIT clustering summary."""
    sizes      = [len(v) for v in clusters.values()]
    singletons = sum(1 for s in sizes if s == 1)
    lines = [
        "=" * 55,
        "  CD-HIT-EST Clustering Report",
        "=" * 55,
        f"  Identity threshold    : {IDENTITY_THRESHOLD * 100:.0f}%",
        f"  Input sequences       : {input_count}",
        f"  Output representatives: {output_count}",
        f"  Sequences removed     : {input_count - output_count}",
        f"  Reduction             : "
        f"{(input_count - output_count) / input_count * 100:.1f}%",
        "",
        "  Cluster statistics:",
        f"    Total clusters      : {len(clusters)}",
        f"    Singleton clusters  : {singletons}",
        f"    Largest cluster     : {max(sizes)} members",
        f"    Average cluster size: {sum(sizes)/len(sizes):.2f}",
        "=" * 55,
    ]
    report_str = '\n'.join(lines)
    print(report_str)
    with open(filename, 'w') as f:
        f.write(report_str + '\n')
    print(f"\n  Report saved to: {filename}")


def report_class_distribution(fasta_file):
    """Reports Class A/B/C/D counts in a FASTA file."""
    counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    for header, _ in read_fasta(fasta_file):
        n = int(header.split('_')[1])
        if   n <= 1250: counts['A'] += 1
        elif n <= 2500: counts['B'] += 1
        elif n <= 3750: counts['C'] += 1
        else:           counts['D'] += 1
    print("\n  Class distribution after CD-HIT:")
    for cls, count in sorted(counts.items()):
        print(f"    Class {cls}: {count}")


if __name__ == "__main__":
    print("=" * 55)
    print("  CD-HIT-EST Redundancy Removal — Script 4")
    print("=" * 55)

    input_count = len(read_fasta(INPUT_FASTA))
    print(f"\n  Input sequences: {input_count}")

    if not run_cdhit(INPUT_FASTA, OUTPUT_PREFIX):
        print("\n  CD-HIT-EST failed. Exiting.")
        exit(1)

    output_count = len(read_fasta(OUTPUT_FASTA))
    clusters, representatives, _ = parse_clstr_file(CLUSTER_FILE)

    write_report(clusters, input_count, output_count)
    report_class_distribution(OUTPUT_FASTA)

    print(f"\n  Non-redundant FASTA : {OUTPUT_FASTA}")
    print(f"  Cluster file        : {CLUSTER_FILE}")




# ─────────────────────────────────────────────────────────────
#   SCRIPT 5 — RNApdist Ensemble Diversity Clustering
# ─────────────────────────────────────────────────────────────
"""
Clusters aptamer sequences by thermodynamic ensemble profile
similarity using RNApdist.

Unlike bp-distance (single MFE structure), RNApdist compares
full base-pair probability matrices from the partition function.
Two sequences may share MFE topology but differ in ensemble —
RNApdist captures this difference.

Distance metric : Euclidean distance between base-pair
                  probability vectors (RNApdist default)
Clustering      : Hierarchical average linkage
Cutoff          : 1.0 — sequences with distance < 1.0 are
                  considered redundant

All pairwise distances are computed in one batch RNApdist
call using -Xm (matrix mode).

Input  : aptamer_cdhit.fasta          (844 sequences from Script 4)
Outputs: aptamer_ensemble_clustered.fasta
         ensemble_clusters.tsv
         ensemble_distance_matrix.npy     ← used by Script 7
         ordered_ids_for_rnadist.txt      ← used by Script 7
         ensemble_distance_distribution.png
"""

import subprocess
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import defaultdict

# ── Parameters ───────────────────────────────────────────────

INPUT_FASTA     = "aptamer_cdhit.fasta"
OUTPUT_FASTA    = "aptamer_ensemble_clustered.fasta"
CLUSTER_REPORT  = "ensemble_clusters.tsv"
DIST_PLOT       = "ensemble_distance_distribution.png"
DISTANCE_CUTOFF = 1.0

# ── Helper functions ─────────────────────────────────────────

def read_fasta(filename):
    """Reads FASTA into ordered list of (header, sequence) tuples."""
    records, header, seq = [], None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.upper())
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def build_rnapdist_matrix(ordered_ids, seq_dict):
    """
    Computes all-vs-all RNApdist distances in a single batch call
    using matrix mode (-Xm). Feeds all sequences via stdin.
    Parses the lower-triangle output into a symmetric numpy matrix.
    Returns numpy array of shape (n, n).
    """
    n           = len(ordered_ids)
    total_pairs = n * (n - 1) // 2
    print(f"  Computing {total_pairs:,} pairwise distances...")
    print(f"  (Single batch call — will complete in minutes)")

    stdin_input = '\n'.join(seq_dict[i] for i in ordered_ids) + '\n'
    result = subprocess.run(
        ['RNApdist', '-Xm'],
        input=stdin_input, capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"RNApdist failed: {result.stderr.strip()}")

    mat   = np.zeros((n, n))
    lines = [l.strip() for l in result.stdout.strip().split('\n')
             if l.strip() and not l.startswith('>')]
    for i, line in enumerate(lines):
        for j, val in enumerate(float(x) for x in line.split()):
            mat[i + 1][j] = val
            mat[j][i + 1] = val

    print(f"  Matrix complete: {mat.shape}")
    return mat


def save_distance_plot(dist_matrix, cutoff, filename):
    """Saves histogram of all pairwise RNApdist distances."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        flat = dist_matrix[np.triu_indices_from(dist_matrix, k=1)]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(flat, bins=50, color='steelblue', edgecolor='white')
        ax.axvline(cutoff, color='red', linestyle='--',
                   linewidth=1.5, label=f'Cutoff = {cutoff}')
        ax.set_xlabel('RNApdist Ensemble Distance', fontsize=12)
        ax.set_ylabel('Pair count', fontsize=12)
        ax.set_title(f'Pairwise RNApdist Distance Distribution '
                     f'({len(dist_matrix)} sequences)', fontsize=13)
        ax.legend()
        fig.tight_layout()
        fig.savefig(filename, dpi=300)
        plt.close()
        print(f"  Distance plot saved: {filename}")
    except ImportError:
        print("  matplotlib not available — skipping distance plot")


def cluster_and_select(dist_matrix, ordered_ids, cutoff):
    """
    Hierarchical average-linkage clustering at given cutoff.
    Representative = member with lowest Aptamer number.
    Returns (representatives list, num_clusters, cluster_members dict).
    """
    if len(ordered_ids) < 2:
        return ordered_ids, 1, {0: ordered_ids}

    condensed = squareform(dist_matrix)
    Z         = linkage(condensed, method='average')
    labels    = fcluster(Z, t=cutoff, criterion='distance')

    cluster_members = defaultdict(list)
    for idx, cid in enumerate(labels):
        cluster_members[cid].append(ordered_ids[idx])

    representatives = []
    for cid in sorted(cluster_members):
        members = sorted(cluster_members[cid],
                         key=lambda x: int(x.split('_')[1]))
        representatives.append(members[0])

    return representatives, len(set(labels)), cluster_members


def report_class_distribution(ids):
    """Reports Class A/B/C/D counts for a list of Aptamer IDs."""
    counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    for apt_id in ids:
        n = int(apt_id.split('_')[1])
        if   n <= 1250: counts['A'] += 1
        elif n <= 2500: counts['B'] += 1
        elif n <= 3750: counts['C'] += 1
        else:           counts['D'] += 1
    for cls, count in sorted(counts.items()):
        print(f"    Class {cls}: {count}")


def write_cluster_report(cluster_members, filename=CLUSTER_REPORT):
    """Writes TSV: cluster_id, representative, member, class."""
    with open(filename, 'w') as f:
        f.write("cluster_id\trepresentative\tmember\tclass\n")
        for cid in sorted(cluster_members):
            members = sorted(cluster_members[cid],
                             key=lambda x: int(x.split('_')[1]))
            rep = members[0]
            for member in members:
                n   = int(member.split('_')[1])
                cls = ('A' if n <= 1250 else 'B' if n <= 2500
                       else 'C' if n <= 3750 else 'D')
                f.write(f"{cid}\t{rep}\t{member}\t{cls}\n")
    print(f"  Cluster report  : {filename}")


if __name__ == "__main__":
    print("=" * 55)
    print("  RNApdist Ensemble Clustering — Script 5")
    print("=" * 55)

    # Step 1 — Load sequences
    print("\nStep 1: Loading sequences...")
    records     = read_fasta(INPUT_FASTA)
    ordered_ids = [r[0] for r in records]
    seq_dict    = {r[0]: r[1] for r in records}
    print(f"  Loaded {len(records)} sequences")

    # Step 2 — Build distance matrix
    print("\nStep 2: Building RNApdist distance matrix...")
    dist_matrix = build_rnapdist_matrix(ordered_ids, seq_dict)
    np.save("ensemble_distance_matrix.npy", dist_matrix)
    print("  Matrix saved    : ensemble_distance_matrix.npy")

    with open("ordered_ids_for_rnadist.txt", 'w') as f:
        f.write('\n'.join(ordered_ids) + '\n')
    print("  Ordered IDs saved: ordered_ids_for_rnadist.txt")

    save_distance_plot(dist_matrix, DISTANCE_CUTOFF, DIST_PLOT)

    # Step 3 — Cluster
    print(f"\nStep 3: Clustering at cutoff = {DISTANCE_CUTOFF}...")
    representatives, num_clusters, cluster_members = cluster_and_select(
        dist_matrix, ordered_ids, DISTANCE_CUTOFF)
    print(f"  Clusters formed      : {num_clusters}")
    print(f"  Representatives kept : {len(representatives)}")

    # Step 4 — Write outputs
    print("\nStep 4: Writing outputs...")
    sorted_reps = sorted(representatives,
                         key=lambda x: int(x.split('_')[1]))
    with open(OUTPUT_FASTA, 'w') as f:
        for apt_id in sorted_reps:
            f.write(f">{apt_id}\n{seq_dict[apt_id]}\n")

    write_cluster_report(cluster_members)

    print(f"\n{'=' * 55}")
    print(f"  Input              : {len(records)} sequences")
    print(f"  Clusters formed    : {num_clusters}")
    print(f"  Representatives    : {len(representatives)}")
    print(f"  Output FASTA       : {OUTPUT_FASTA}")
    print(f"\n  Class distribution after ensemble clustering:")
    report_class_distribution(representatives)
    print(f"{'=' * 55}")




# ─────────────────────────────────────────────────────────────
#   SCRIPT 6 — Final Library Preparation
# ─────────────────────────────────────────────────────────────
"""
Compiles the final 257 aptamer candidates into a clean,
annotated library ready for 3D structure prediction and docking.

Pulls together all metrics from Scripts 2–5:
    sequence, length, class, MFE total, MFE per nt,
    ensemble diversity, loop length, GGG count, stem pairs,
    ensemble cluster ID.

Inputs  : aptamer_ensemble_clustered.fasta
          filter_report.tsv
          ensemble_clusters.tsv
Outputs : final_aptamer_library.fasta
          final_aptamer_metadata.tsv
          final_library_summary.txt
"""

import re
from collections import defaultdict

# ── Parameters ───────────────────────────────────────────────

INPUT_FASTA     = "aptamer_ensemble_clustered.fasta"
FILTER_REPORT   = "filter_report.tsv"
CLUSTER_REPORT  = "ensemble_clusters.tsv"
OUTPUT_FASTA    = "final_aptamer_library.fasta"
OUTPUT_METADATA = "final_aptamer_metadata.tsv"
OUTPUT_SUMMARY  = "final_library_summary.txt"

# ── Helper functions ─────────────────────────────────────────

def read_fasta(filename):
    """Reads FASTA into ordered list of (header, sequence) tuples."""
    records, header, seq = [], None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.upper())
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def get_class(aptamer_id):
    """Determines design class from Aptamer number (1–5000)."""
    n = int(aptamer_id.split('_')[1])
    if   n <= 1250: return 'A'
    elif n <= 2500: return 'B'
    elif n <= 3750: return 'C'
    else:           return 'D'


def class_description(cls):
    return {'A': 'G4-Core loop', 'B': 'UGGG-Anchor loop',
            'C': 'Distributed-G loop', 'D': 'Mixed scaffold loop'}.get(cls, 'Unknown')


def load_filter_report(filename):
    """Loads filter_report.tsv → dict: aptamer_id → metrics dict."""
    metrics = {}
    with open(filename, 'r') as f:
        headers = f.readline().strip().split('\t')
        for line in f:
            row = dict(zip(headers, line.strip().split('\t')))
            metrics[row['aptamer_id']] = row
    return metrics


def load_cluster_report(filename):
    """Loads ensemble_clusters.tsv → dict: aptamer_id → cluster_id."""
    membership = {}
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                membership[parts[2]] = parts[0]
    return membership


# ── Main compilation ─────────────────────────────────────────

def compile_library(records, filter_metrics, cluster_membership):
    """
    Assembles one dict per aptamer containing all metrics.
    Returns list of dicts sorted by original Aptamer number.
    """
    compiled = []
    for original_id, seq in records:
        cls = get_class(original_id)
        fm  = filter_metrics.get(original_id, {})
        compiled.append({
            'original_id' : original_id,
            'class'       : cls,
            'class_desc'  : class_description(cls),
            'sequence'    : seq,
            'length'      : len(seq),
            'mfe_total'   : fm.get('mfe_total',   'NA'),
            'mfe_per_nt'  : fm.get('mfe_per_nt',  'NA'),
            'diversity'   : fm.get('diversity',   'NA'),
            'loop_len'    : fm.get('loop_len',    'NA'),
            'stem_pairs'  : fm.get('stem_pairs',  'NA'),
            'ggg_in_loop' : fm.get('ggg_in_loop', 'NA'),
            'cluster_id'  : cluster_membership.get(original_id, 'NA'),
        })
    return compiled


# ── Output writers ───────────────────────────────────────────

def write_final_fasta(compiled, filename=OUTPUT_FASTA):
    with open(filename, 'w') as f:
        for e in compiled:
            f.write(f">{e['original_id']}\n{e['sequence']}\n")
    print(f"  Final FASTA    : {filename}  ({len(compiled)} sequences)")


def write_metadata(compiled, filename=OUTPUT_METADATA):
    """Writes full metadata TSV — becomes Table S1 in supplementary."""
    cols = ['original_id', 'class', 'class_desc', 'sequence', 'length',
            'mfe_total', 'mfe_per_nt', 'diversity', 'loop_len',
            'stem_pairs', 'ggg_in_loop', 'cluster_id']
    with open(filename, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for e in compiled:
            f.write('\t'.join(str(e.get(c, 'NA')) for c in cols) + '\n')
    print(f"  Metadata TSV   : {filename}")


def write_summary(compiled, filename=OUTPUT_SUMMARY):
    """
    Writes a human-readable final library summary.
    This text maps directly to the Methods section paragraph
    describing the computational library preparation.
    """
    total         = len(compiled)
    class_counts  = defaultdict(int)
    length_counts = defaultdict(int)
    mfe_vals, div_vals = [], []

    for e in compiled:
        class_counts[e['class']] += 1
        length_counts[e['length']] += 1
        try:
            mfe_vals.append(float(e['mfe_total']))
            div_vals.append(float(e['diversity']))
        except (ValueError, TypeError):
            pass

    lines = [
        "=" * 60,
        "  FINAL APTAMER LIBRARY — SUMMARY REPORT",
        "=" * 60,
        "",
        "  LIBRARY COMPOSITION",
        f"    Total aptamers          : {total}",
        f"    Class A (G4-Core)       : {class_counts['A']}",
        f"    Class B (UGGG-Anchor)   : {class_counts['B']}",
        f"    Class C (Distributed-G) : {class_counts['C']}",
        f"    Class D (Mixed scaffold): {class_counts['D']}",
        "",
        "  LENGTH DISTRIBUTION",
    ]
    for ln in sorted(length_counts):
        lines.append(f"    {ln} nt : {length_counts[ln]} sequences")

    lines += [
        "",
        "  BIOPHYSICAL METRICS",
        f"    MFE range    : {min(mfe_vals):.2f} to {max(mfe_vals):.2f} kcal/mol",
        f"    MFE mean     : {sum(mfe_vals)/len(mfe_vals):.2f} kcal/mol",
        f"    Diversity mean: {sum(div_vals)/len(div_vals):.2f}",
        "",
        "  PIPELINE FUNNEL",
        "    Script 1 — Generated          : 5,000",
        "    Script 2 — RNAfold folded      : 4,923",
        "    Script 3 — Biophysical filter  : 1,179",
        "    Script 4 — CD-HIT 80%          :   844",
        "    Script 5 — RNApdist cutoff 1.0 :   257",
        f"    Script 6 — Final library       :   {total}",
        "",
        "  TARGET PHOSPHORYLATION SITES",
        "    Early-stage  : pT181, pT217, pT231, pS202",
        "    MT-binding   : pS262, pS305, pS356",
        "    Aggregation  : pS396, pS404, pS409, pS412, pS422",
        "",
        "=" * 60,
        "  Ready for 3D structure prediction and docking.",
        "=" * 60,
    ]
    report_str = '\n'.join(lines)
    print(report_str)
    with open(filename, 'w') as f:
        f.write(report_str + '\n')
    print(f"\n  Summary        : {filename}")


if __name__ == "__main__":
    print("=" * 60)
    print("  Final Library Preparation — Script 6")
    print("=" * 60)

    records            = read_fasta(INPUT_FASTA)
    filter_metrics     = load_filter_report(FILTER_REPORT)
    cluster_membership = load_cluster_report(CLUSTER_REPORT)

    print(f"\n  Sequences loaded  : {len(records)}")
    print(f"  Filter metrics    : {len(filter_metrics)}")
    print(f"  Cluster data      : {len(cluster_membership)}")

    compiled = compile_library(records, filter_metrics, cluster_membership)

    write_final_fasta(compiled)
    write_metadata(compiled)
    write_summary(compiled)

    print("\n" + "=" * 60)
    print("  Script 6 Complete.")
    print("  final_aptamer_library.fasta ready for")
    print("  3D structure prediction (trRosettaRNA2).")
    print("=" * 60)




# ─────────────────────────────────────────────────────────────
#   SCRIPT 7 — Figures and Visualization
# ─────────────────────────────────────────────────────────────
"""
Generates all publication-quality figures for the aptamer
design pipeline thesis chapter.

Figures produced:
    Fig 01 — Pipeline Funnel (sequence attrition)
    Fig 02 — Class Distribution Across Pipeline Stages
    Fig 03 — MFE Total Distribution (final 257)
    Fig 04 — MFE per Nucleotide Distribution
    Fig 05 — Ensemble Diversity Distribution
    Fig 06 — MFE vs Diversity Scatter (colored by class)
    Fig 07 — MFE per nt vs Stem Pairs Scatter
    Fig 08 — Loop Length Distribution by Class
    Fig 09 — Stem Pairs Distribution by Class
    Fig 10 — GGG-in-Loop Distribution by Class
    Fig 11 — Aptamer Length Distribution by Class
    Fig 12 — RNApdist Pairwise Distance Heatmap (257×257)
    Fig 13 — Ensemble Cluster Size Distribution
    Fig 14 — Hierarchical Clustering Dendrogram (844 aptamers)
    Fig 15 — Pairwise RNApdist Distance Distribution

Inputs:
    filter_report.tsv                (Script 3)
    final_aptamer_metadata.tsv       (Script 6)
    ensemble_clusters.tsv            (Script 5)
    aptamer_ensemble_clustered.fasta (Script 5)
    ensemble_distance_matrix.npy     (Script 5)
    ordered_ids_for_rnadist.txt      (Script 5)

Outputs:
    fig01_pipeline_funnel.png
    fig02_class_distribution_stages.png
    fig03_mfe_total_distribution.png
    fig04_mfe_per_nt_distribution.png
    fig05_diversity_distribution.png
    fig06_mfe_vs_diversity_scatter.png
    fig07_mfe_per_nt_vs_stem_pairs.png
    fig08_loop_length_distribution.png
    fig09_stem_pairs_distribution.png
    fig10_ggg_in_loop_distribution.png
    fig11_length_distribution.png
    fig12_rnapdist_heatmap.png
    fig13_cluster_size_distribution.png
    fig14_dendrogram.png
    fig15_distance_distribution.png
"""

import os
import csv
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# ── Global settings ──────────────────────────────────────────

DPI            = 300
FIGSIZE_SINGLE = (10, 6)
FIGSIZE_WIDE   = (14, 6)
FIGSIZE_SQUARE = (8, 8)

CLASS_COLORS = {
    'A': '#2E86AB',   # Steel blue   — G4-Core
    'B': '#A23B72',   # Purple-rose  — UGGG-Anchor
    'C': '#F18F01',   # Amber        — Distributed-G
    'D': '#C73E1D',   # Rust red     — Mixed scaffold
}
CLASS_LABELS = {
    'A': 'Class A — G4-Core',
    'B': 'Class B — UGGG-Anchor',
    'C': 'Class C — Distributed-G',
    'D': 'Class D — Mixed scaffold',
}

# ── Input files ──────────────────────────────────────────────

FILTER_REPORT  = "filter_report.tsv"
METADATA_TSV   = "final_aptamer_metadata.tsv"
CLUSTER_REPORT = "ensemble_clusters.tsv"
FINAL_FASTA    = "aptamer_ensemble_clustered.fasta"
MATRIX_FILE    = "ensemble_distance_matrix.npy"
ORDERED_IDS    = "ordered_ids_for_rnadist.txt"

# ── Data loaders ─────────────────────────────────────────────

def load_metadata(filename):
    """Loads final_aptamer_metadata.tsv into a list of dicts."""
    with open(filename, 'r') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_filter_report(filename, keep_ids):
    """
    Loads filter_report.tsv, keeping only the final 257 aptamers.
    Attaches 'class' key to each row.
    """
    records = []
    with open(filename, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['aptamer_id'] in keep_ids:
                row['class'] = get_class(row['aptamer_id'])
                records.append(row)
    return records


def load_cluster_report(filename):
    """Loads ensemble_clusters.tsv → dict: cluster_id → [members]."""
    clusters = defaultdict(list)
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                clusters[parts[0]].append(parts[2])
    return clusters


def get_class(aptamer_id):
    n = int(aptamer_id.split('_')[1])
    if   n <= 1250: return 'A'
    elif n <= 2500: return 'B'
    elif n <= 3750: return 'C'
    else:           return 'D'


def safe_float(val, default=None):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


# ── Figure 01 — Pipeline Funnel ──────────────────────────────

def fig01_pipeline_funnel():
    """Horizontal funnel showing sequence attrition at each step."""
    steps  = ['Script 1\nGenerated', 'Script 2\nRNAfold Filter',
               'Script 3\nBiophysical Filter', 'Script 4\nCD-HIT 80%',
               'Script 5\nRNApdist Cluster']
    counts = [5000, 4923, 1179, 844, 257]
    colors = ['#1a535c', '#4ecdc4', '#ff6b6b', '#ffe66d', '#2E86AB']

    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    max_val = max(counts)

    for i, (step, count, color) in enumerate(
            zip(reversed(steps), reversed(counts), reversed(colors))):
        w = count / max_val
        l = (1 - w) / 2
        ax.barh(i, w, left=l, height=0.5,
                color=color, edgecolor='white', linewidth=1.5)
        ax.text(0.5, i, f'{count:,}', ha='center', va='center',
                fontsize=14, fontweight='bold', color='white')
        ax.text(l - 0.01, i, step, ha='right', va='center',
                fontsize=11, color='#333333')

    for i, (count, prev) in enumerate(
            zip(reversed(counts), [None] + list(reversed(counts))[:-1])):
        if prev is not None:
            ax.text(0.98, i, f'↓ {count / prev * 100:.1f}%',
                    ha='right', va='center', fontsize=10,
                    color='#666666', style='italic')

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(steps) - 0.5)
    ax.axis('off')
    ax.set_title('Aptamer Design Pipeline — Sequence Attrition',
                 fontsize=15, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('fig01_pipeline_funnel.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 01 saved.")


# ── Figure 02 — Class Distribution Across Pipeline Stages ────

def fig02_class_distribution_stages():
    """Grouped bar chart of Class A/B/C/D counts at each stage."""
    stages = ['Generated\n(Script 1)', 'Biophysical\n(Script 3)',
              'CD-HIT\n(Script 4)', 'RNApdist\n(Script 5)',
              'Final\n(Script 6)']
    data   = {
        'A': [1250, 619, 457, 127, 127],
        'B': [1250, 303, 218,  61,  61],
        'C': [1250, 147, 108,  39,  39],
        'D': [1250, 110,  61,  30,  30],
    }
    x, width = np.arange(len(stages)), 0.18
    fig, ax  = plt.subplots(figsize=FIGSIZE_WIDE)

    for i, cls in enumerate(['A', 'B', 'C', 'D']):
        bars = ax.bar(x + (i - 1.5) * width, data[cls], width,
                      label=CLASS_LABELS[cls],
                      color=CLASS_COLORS[cls], edgecolor='white')
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2,
                    h + 8, str(int(h)), ha='center', va='bottom',
                    fontsize=7.5, color='#333333')

    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=11)
    ax.set_ylabel('Number of Aptamers', fontsize=12)
    ax.set_title('Class Distribution Across Pipeline Stages',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(0, 1500)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig('fig02_class_distribution_stages.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 02 saved.")


# ── Figures 03–05 — MFE / Diversity Distributions ────────────

def _histogram_by_class(records, key, xlabel, title, filename):
    """Generic helper for per-class overlapping histograms."""
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    for cls in ['A', 'B', 'C', 'D']:
        vals = [safe_float(r[key]) for r in records
                if r['class'] == cls and safe_float(r[key]) is not None]
        ax.hist(vals, bins=20, alpha=0.65, color=CLASS_COLORS[cls],
                label=CLASS_LABELS[cls], edgecolor='white')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  {filename} saved.")


def fig03_mfe_total_distribution(records):
    _histogram_by_class(records, 'mfe_total', 'MFE (kcal/mol)',
                        'MFE Distribution — Final 257 Aptamers',
                        'fig03_mfe_total_distribution.png')


def fig04_mfe_per_nt_distribution(records):
    _histogram_by_class(records, 'mfe_per_nt', 'MFE per Nucleotide (kcal/mol/nt)',
                        'MFE per Nucleotide — Final 257 Aptamers',
                        'fig04_mfe_per_nt_distribution.png')


def fig05_diversity_distribution(records):
    _histogram_by_class(records, 'diversity', 'Ensemble Diversity',
                        'Ensemble Diversity — Final 257 Aptamers',
                        'fig05_diversity_distribution.png')


# ── Figures 06–07 — Scatter Plots ────────────────────────────

def _scatter_by_class(records, x_key, y_key, xlabel, ylabel, title, filename):
    """Generic helper for per-class scatter plots."""
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    for cls in ['A', 'B', 'C', 'D']:
        xs = [safe_float(r[x_key]) for r in records if r['class'] == cls
              and safe_float(r[x_key]) is not None and safe_float(r[y_key]) is not None]
        ys = [safe_float(r[y_key]) for r in records if r['class'] == cls
              and safe_float(r[x_key]) is not None and safe_float(r[y_key]) is not None]
        ax.scatter(xs, ys, alpha=0.7, s=40, color=CLASS_COLORS[cls],
                   label=CLASS_LABELS[cls], edgecolors='white', linewidths=0.4)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3)
    ax.xaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  {filename} saved.")


def fig06_mfe_vs_diversity_scatter(records):
    _scatter_by_class(records, 'mfe_total', 'diversity',
                      'MFE (kcal/mol)', 'Ensemble Diversity',
                      'Thermodynamic Landscape — Final 257 Aptamers',
                      'fig06_mfe_vs_diversity_scatter.png')


def fig07_mfe_per_nt_vs_stem_pairs(records):
    _scatter_by_class(records, 'stem_pairs', 'mfe_per_nt',
                      'Stem Pairs (count)', 'MFE per nt (kcal/mol/nt)',
                      'Structural Stability vs Stem Pairs',
                      'fig07_mfe_per_nt_vs_stem_pairs.png')


# ── Figures 08–11 — Structural Feature Distributions ─────────

def _grouped_bar_by_class(records, key, xlabel, title, filename):
    """Generic grouped bar for discrete structural features."""
    all_vals = sorted(set(int(safe_float(r[key])) for r in records
                          if safe_float(r[key]) is not None))
    x, width = np.arange(len(all_vals)), 0.18
    fig, ax  = plt.subplots(figsize=FIGSIZE_SINGLE)

    for i, cls in enumerate(['A', 'B', 'C', 'D']):
        counts = [sum(1 for r in records if r['class'] == cls
                      and safe_float(r[key]) == v) for v in all_vals]
        ax.bar(x + (i - 1.5) * width, counts, width,
               label=CLASS_LABELS[cls],
               color=CLASS_COLORS[cls], edgecolor='white')

    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in all_vals])
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  {filename} saved.")


def fig08_loop_length_distribution(records):
    _grouped_bar_by_class(records, 'loop_len', 'Loop Length (nt)',
                          'Loop Length Distribution by Class',
                          'fig08_loop_length_distribution.png')


def fig09_stem_pairs_distribution(records):
    _grouped_bar_by_class(records, 'stem_pairs', 'Stem Pairs (count)',
                          'Stem Pair Distribution by Class',
                          'fig09_stem_pairs_distribution.png')


def fig10_ggg_in_loop_distribution(records):
    _grouped_bar_by_class(records, 'ggg_in_loop', 'GGG Count in Loop',
                          'GGG Motif Frequency in Loop by Class',
                          'fig10_ggg_in_loop_distribution.png')


def fig11_length_distribution(records):
    """Length distribution (28/32/36/40 nt) grouped by class."""
    lengths  = [28, 32, 36, 40]
    x, width = np.arange(len(lengths)), 0.18
    fig, ax  = plt.subplots(figsize=FIGSIZE_SINGLE)

    for i, cls in enumerate(['A', 'B', 'C', 'D']):
        counts = [sum(1 for r in records if r['class'] == cls
                      and safe_float(r['length']) == ln) for ln in lengths]
        bars = ax.bar(x + (i - 1.5) * width, counts, width,
                      label=CLASS_LABELS[cls],
                      color=CLASS_COLORS[cls], edgecolor='white')
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        h + 0.3, str(int(h)),
                        ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([f'{ln} nt' for ln in lengths], fontsize=12)
    ax.set_xlabel('Aptamer Length (nt)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Aptamer Length Distribution by Class',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig('fig11_length_distribution.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 11 saved.")


# ── Figure 12 — RNApdist Heatmap (257×257) ───────────────────

def fig12_rnapdist_heatmap():
    """
    Loads pre-saved distance matrix from Script 5 and plots
    a 257×257 pairwise heatmap sorted by class.
    Requires: ensemble_distance_matrix.npy, ordered_ids_for_rnadist.txt
    """
    if not os.path.exists(MATRIX_FILE):
        print(f"  Fig 12 SKIPPED: {MATRIX_FILE} not found.")
        return

    dist_matrix = np.load(MATRIX_FILE)
    final_ids   = [l.strip()[1:] for l in open(FINAL_FASTA)
                   if l.startswith('>')]

    with open(ORDERED_IDS, 'r') as f:
        all_ids = [l.strip() for l in f if l.strip()]

    id_to_idx        = {aid: i for i, aid in enumerate(all_ids)}
    final_idx        = [id_to_idx[aid] for aid in final_ids if aid in id_to_idx]
    final_ids_matched = [all_ids[i] for i in final_idx]
    sub_matrix       = dist_matrix[np.ix_(final_idx, final_idx)]

    n          = len(final_ids_matched)
    classes    = [get_class(aid) for aid in final_ids_matched]
    sort_order = sorted(range(n), key=lambda i: classes[i])
    mat_sorted = sub_matrix[np.ix_(sort_order, sort_order)]
    cls_sorted = [classes[i] for i in sort_order]

    cmap = LinearSegmentedColormap.from_list(
        'rnapdist', ['#ffffff', '#4ecdc4', '#1a535c'])
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(mat_sorted, cmap=cmap, aspect='auto',
                   interpolation='nearest')
    plt.colorbar(im, ax=ax, label='RNApdist Ensemble Distance', shrink=0.8)

    boundaries, prev_cls = {}, None
    for i, cls in enumerate(cls_sorted):
        if cls != prev_cls:
            if i > 0:
                ax.axhline(i - 0.5, color='red', linewidth=1.2, alpha=0.8)
                ax.axvline(i - 0.5, color='red', linewidth=1.2, alpha=0.8)
            boundaries[cls] = i
            prev_cls = cls
    
    class_order = list(boundaries.keys())
    for k, cls in enumerate(class_order):
        start = boundaries[cls]
        end   = boundaries[class_order[k + 1]] if k + 1 < len(class_order) else n
        ax.text((start + end) / 2, n + 6, f'Class {cls}',
                ha='center', va='top', fontsize=10,
                color=CLASS_COLORS[cls], fontweight='bold')

    ax.set_title(
        'RNApdist Pairwise Distance Heatmap — Final 257 Aptamers\n'
        '(sorted by class; computed from base-pair probability profiles)',
        fontsize=13, fontweight='bold')
    ax.set_xlabel('Aptamer Index (sorted by class)', fontsize=11)
    ax.set_ylabel('Aptamer Index (sorted by class)', fontsize=11)
    plt.tight_layout()
    plt.savefig('fig12_rnapdist_heatmap.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 12 saved.")


# ── Figure 13 — Cluster Size Distribution ────────────────────

def fig13_cluster_size_distribution(cluster_data):
    """Bar chart + histogram of ensemble cluster sizes from Script 5."""
    sizes = sorted([len(v) for v in cluster_data.values()], reverse=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGSIZE_WIDE)

    ax1.bar(range(1, len(sizes) + 1), sizes, color='#2E86AB',
            edgecolor='white', linewidth=0.5)
    ax1.set_xlabel('Cluster Index (sorted by size)', fontsize=11)
    ax1.set_ylabel('Sequences in Cluster', fontsize=11)
    ax1.set_title('Cluster Sizes (all 257 clusters)',
                  fontsize=12, fontweight='bold')
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.yaxis.grid(True, alpha=0.3)
    ax1.set_axisbelow(True)

    ax2.hist(sizes, bins=max(sizes), color='#A23B72',
             edgecolor='white', linewidth=0.5)
    ax2.set_xlabel('Cluster Size', fontsize=11)
    ax2.set_ylabel('Number of Clusters', fontsize=11)
    ax2.set_title('Cluster Size Frequency', fontsize=12, fontweight='bold')
    ax2.spines[['top', 'right']].set_visible(False)
    ax2.yaxis.grid(True, alpha=0.3)
    ax2.set_axisbelow(True)
    ax2.text(0.97, 0.95,
             f'Singletons: {sum(1 for s in sizes if s == 1)}\n'
             f'Max size: {max(sizes)}\nMedian: {int(np.median(sizes))}',
             transform=ax2.transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.4))

    plt.suptitle('Ensemble Cluster Analysis — Script 5 Output',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig('fig13_cluster_size_distribution.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 13 saved.")


# ── Figure 14 — Hierarchical Clustering Dendrogram ───────────

def fig14_dendrogram():
    """
    Rebuilds the hierarchical clustering dendrogram from the saved
    distance matrix. Leaf colors encode class membership.
    Requires: ensemble_distance_matrix.npy, ordered_ids_for_rnadist.txt
    """
    if not os.path.exists(MATRIX_FILE):
        print("  Fig 14 SKIPPED: ensemble_distance_matrix.npy not found.")
        return

    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

    CUTOFF      = 1.0
    dist_matrix = np.load(MATRIX_FILE)

    with open(ORDERED_IDS, 'r') as f:
        all_ids = [l.strip() for l in f if l.strip()]

    condensed = squareform(dist_matrix)
    Z         = linkage(condensed, method='average')

    fig, ax = plt.subplots(figsize=(20, 8))
    dendrogram(Z, no_labels=True, color_threshold=CUTOFF,
               above_threshold_color='#AAAAAA', ax=ax)

    ax.axhline(y=CUTOFF, color='red', linestyle='--',
               linewidth=1.5, label=f'Cutoff = {CUTOFF}')
    ax.set_title(
        'Hierarchical Clustering Dendrogram — RNApdist Ensemble Distances\n'
        '(844 aptamers after CD-HIT; cutoff = 1.0 → 257 clusters)',
        fontsize=14, fontweight='bold')
    ax.set_xlabel('Aptamer Candidates', fontsize=12)
    ax.set_ylabel('RNApdist Ensemble Distance', fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)

    patches = [mpatches.Patch(color=CLASS_COLORS[c], label=CLASS_LABELS[c])
               for c in ['A', 'B', 'C', 'D']]
    patches.append(mpatches.Patch(color='red', label=f'Cutoff = {CUTOFF}'))
    ax.legend(handles=patches, fontsize=10, loc='upper right')

    plt.tight_layout()
    plt.savefig('fig14_dendrogram.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 14 saved.")


# ── Figure 15 — Pairwise Distance Distribution ───────────────

def fig15_distance_distribution():
    """
    Histogram of all 355,746 pairwise RNApdist distances from the
    844-sequence matrix. Shows where cutoff = 1.0 sits.
    Requires: ensemble_distance_matrix.npy
    """
    if not os.path.exists(MATRIX_FILE):
        print("  Fig 15 SKIPPED: ensemble_distance_matrix.npy not found.")
        return

    CUTOFF      = 1.0
    dist_matrix = np.load(MATRIX_FILE)
    n           = dist_matrix.shape[0]
    upper_tri   = dist_matrix[np.triu_indices(n, k=1)]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.hist(upper_tri, bins=80, color='#2E86AB',
            edgecolor='white', linewidth=0.4, alpha=0.85)
    ax.axvline(x=CUTOFF, color='red', linestyle='--',
               linewidth=2.0, label=f'Cutoff = {CUTOFF}')
    ax.axvspan(0, CUTOFF, alpha=0.12, color='red',
               label='Clustered together (< 1.0)')
    ax.set_xlabel('RNApdist Ensemble Distance', fontsize=12)
    ax.set_ylabel('Pair Count', fontsize=12)
    ax.set_title(
        'Pairwise RNApdist Distance Distribution\n'
        '(355,746 pairs from 844 aptamers)',
        fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig('fig15_distance_distribution.png', dpi=DPI,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Fig 15 saved.")


# ── Entry point ──────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("  Aptamer Figures — Script 7")
    print("=" * 60)

    print("\nLoading data...")
    metadata       = load_metadata(METADATA_TSV)
    keep_ids       = {r['original_id'] for r in metadata}
    filter_records = load_filter_report(FILTER_REPORT, keep_ids)
    cluster_data   = load_cluster_report(CLUSTER_REPORT)

    print(f"  Metadata       : {len(metadata)} records")
    print(f"  Filter records : {len(filter_records)}")
    print(f"  Clusters       : {len(cluster_data)}")

    print("\nGenerating figures...")
    fig01_pipeline_funnel()
    fig02_class_distribution_stages()
    fig03_mfe_total_distribution(filter_records)
    fig04_mfe_per_nt_distribution(filter_records)
    fig05_diversity_distribution(filter_records)
    fig06_mfe_vs_diversity_scatter(filter_records)
    fig07_mfe_per_nt_vs_stem_pairs(filter_records)
    fig08_loop_length_distribution(filter_records)
    fig09_stem_pairs_distribution(filter_records)
    fig10_ggg_in_loop_distribution(filter_records)
    fig11_length_distribution(filter_records)
    fig12_rnapdist_heatmap()
    fig13_cluster_size_distribution(cluster_data)
    fig14_dendrogram()
    fig15_distance_distribution()

    print("\n" + "=" * 60)
    print("  Script 7 Complete — 15 figures generated.")
    print("=" * 60)
