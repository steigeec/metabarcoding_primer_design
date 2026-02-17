#!/usr/bin/env python3
"""
design_primers.py

Second stage of the anuran metabarcoding pipeline.

Takes large FASTA files produced by fetch_anuran_loci.py and outputs candidate
degenerate primer pairs per locus. Constraint: total amplicon ≤ 65 bp
(including primers) for short-read/degraded-DNA sequencing.

Requires:
  mafft        (conda install -c bioconda mafft)
  trimal       (conda install -c bioconda trimal)
  DEGEPRIME    (git clone https://github.com/envgen/DEGEPRIME)
  biopython    (pip install biopython)

Usage:
  python design_primers.py --degeprime ./DEGEPRIME/DegePrime.pl
  python design_primers.py --degeprime ./DEGEPRIME/DegePrime.pl --loci 12S --seqs-per-family 2
"""

import os
import re
import sys
import csv
import shutil
import random
import bisect
import argparse
import subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqUtils.MeltingTemp import Tm_NN

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALL_LOCI = ["12S", "16S", "cytb", "COI"]

# IUPAC complement table (uppercase; handles ambiguity codes)
_IUPAC_COMP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N', '-': '-',
}

# IUPAC ambiguity → set of bases (used for 3′-end dimer checks)
_IUPAC_BASES = {
    'A': {'A'}, 'T': {'T'}, 'G': {'G'}, 'C': {'C'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'G', 'C'}, '-': set(),
}

# ---------------------------------------------------------------------------
# Genus extraction — copied verbatim from fetch_anuran_loci.py
# ---------------------------------------------------------------------------

_GENUS_RE = re.compile(r'(?:^|\s)([A-Z][a-z]{2,})(?=\s+[a-z])')


def _extract_genus(description: str) -> str:
    """
    Extract the genus name from an NCBI FASTA description line.
    rec.description includes the accession as the first token, so we scan
    for the first word that looks like a genus (Capitalised, ≥3 chars,
    followed by a lowercase word).
    """
    match = _GENUS_RE.search(description)
    return match.group(1) if match else "Unknown"


# ---------------------------------------------------------------------------
# 1. IUPAC utilities
# ---------------------------------------------------------------------------

def iupac_complement(seq: str) -> str:
    """Return the IUPAC complement of *seq* (does NOT reverse)."""
    return ''.join(_IUPAC_COMP.get(b.upper(), 'N') for b in seq)


def iupac_to_bases(char: str) -> set:
    """Return the set of unambiguous bases represented by an IUPAC character."""
    return _IUPAC_BASES.get(char.upper(), {'N'})


# ---------------------------------------------------------------------------
# 2. Validate degeneracy
# ---------------------------------------------------------------------------

def validate_degeneracy(d: int) -> None:
    """
    DEGEPRIME requires degeneracy to be of the form 2^i × 3^j.
    Raises ValueError with suggestions if *d* does not satisfy this.
    """
    n = d
    for p in (2, 3):
        while n % p == 0:
            n //= p
    if n != 1:
        # Generate nearby valid values for the error message
        valid = sorted(
            {2**i * 3**j for i in range(8) for j in range(6)
             if 1 <= 2**i * 3**j <= 1000}
        )
        below = [v for v in valid if v < d]
        above = [v for v in valid if v > d]
        suggestion = []
        if below:
            suggestion.append(str(below[-1]))
        if above:
            suggestion.append(str(above[0]))
        raise ValueError(
            f"--max-degeneracy {d} is not of the form 2^i × 3^j. "
            f"Try: {' or '.join(suggestion)}"
        )


# ---------------------------------------------------------------------------
# 3. Check dependencies
# ---------------------------------------------------------------------------

def check_dependencies(degeprime_path: str) -> None:
    """
    Verify that mafft, trimal, perl, degeprime file, and biopython are available.
    Prints actionable install hints and exits on failure.
    """
    errors = []

    if shutil.which("mafft") is None:
        errors.append(
            "mafft not found. Install with:\n"
            "  conda install -c bioconda mafft"
        )
    if shutil.which("trimal") is None:
        errors.append(
            "trimal not found. Install with:\n"
            "  conda install -c bioconda trimal"
        )
    if shutil.which("perl") is None:
        errors.append("perl not found. Install Perl 5.")

    if not os.path.isfile(degeprime_path):
        errors.append(
            f"DegePrime.pl not found at: {degeprime_path}\n"
            "  Clone with: git clone https://github.com/envgen/DEGEPRIME"
        )

    try:
        from Bio.SeqUtils.MeltingTemp import Tm_NN  # noqa: F401
    except ImportError:
        errors.append(
            "biopython not found. Install with:\n"
            "  pip install biopython"
        )

    if errors:
        print("\nDependency check failed:\n")
        for e in errors:
            print(f"  ERROR: {e}\n")
        sys.exit(1)


# ---------------------------------------------------------------------------
# 4. Parse coverage file
# ---------------------------------------------------------------------------

def parse_coverage_file(path: str) -> dict:
    """
    Parse a *_tax_coverage.txt produced by fetch_anuran_loci.py.
    Returns {family: count} skipping header/dashes lines.
    The count is the last whitespace-delimited token on each data line.
    """
    family_counts = {}
    if not os.path.exists(path):
        return family_counts
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith('=') or line.startswith('-'):
                continue
            tokens = line.split()
            if len(tokens) < 2:
                continue
            # Last token must be a number
            try:
                count = int(tokens[-1].replace(',', ''))
            except ValueError:
                continue
            family = ' '.join(tokens[:-1])
            # Skip header-like lines
            if family.lower() in ('family', 'total sequences', 'families found'):
                continue
            family_counts[family] = count
    return family_counts


# ---------------------------------------------------------------------------
# 5. Subsample records
# ---------------------------------------------------------------------------

def subsample_records(
    fasta_path: str,
    coverage_path: str,
    seqs_per_family: int,
    seed: int,
) -> list:
    """
    Load FASTA, drop sequences < 50 bp or all-N, group by genus, sample up
    to *seqs_per_family* per genus.  Returns a list of SeqRecord objects.
    """
    rng = random.Random(seed)

    print(f"  Loading sequences from {fasta_path} ...")
    all_records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"  {len(all_records):,} records loaded.")

    # Filter short / all-N
    filtered = []
    for rec in all_records:
        seq_str = str(rec.seq).upper()
        if len(seq_str) < 50:
            continue
        if all(b in ('N', '-', '?') for b in seq_str):
            continue
        filtered.append(rec)
    print(f"  {len(filtered):,} records after quality filter (≥50 bp, not all-N).")

    # Group by genus
    by_genus: dict = defaultdict(list)
    for rec in filtered:
        genus = _extract_genus(rec.description)
        by_genus[genus].append(rec)

    n_genera = len(by_genus)
    print(f"  {n_genera:,} unique genera found.")

    # Subsample
    sampled = []
    for genus, recs in sorted(by_genus.items()):
        if len(recs) <= seqs_per_family:
            sampled.extend(recs)
        else:
            sampled.extend(rng.sample(recs, seqs_per_family))

    print(f"  {len(sampled):,} sequences selected (≤{seqs_per_family} per genus).")

    # Log family coverage for reference (informational only)
    family_counts = parse_coverage_file(coverage_path)
    if family_counts:
        print(f"  Reference coverage file lists {len(family_counts):,} families.")

    return sampled


# ---------------------------------------------------------------------------
# 6. Alignment wrappers
# ---------------------------------------------------------------------------

def run_mafft(in_path: str, out_path: str, threads: int) -> None:
    """Align sequences with MAFFT --auto."""
    cmd = ["mafft", "--auto", "--thread", str(threads), in_path]
    print(f"  Running: {' '.join(cmd)}")
    with open(out_path, "w") as fh:
        result = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed (exit {result.returncode}):\n{result.stderr[-2000:]}"
        )


def run_trimal(in_path: str, out_path: str, gt: float) -> None:
    """Trim alignment columns with TrimAl."""
    cmd = ["trimal", "-in", in_path, "-out", out_path, f"-gt", str(gt), "-fasta"]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"TrimAl failed (exit {result.returncode}):\n{result.stderr[-2000:]}"
        )


# ---------------------------------------------------------------------------
# 7. Reverse-complement alignment
# ---------------------------------------------------------------------------

def reverse_complement_alignment(in_path: str, out_path: str) -> None:
    """Write reverse complement of every record in the trimmed FASTA."""
    records = list(SeqIO.parse(in_path, "fasta"))
    rc_records = []
    for rec in records:
        rc = rec.reverse_complement(id=rec.id, description=rec.description)
        rc_records.append(rc)
    SeqIO.write(rc_records, out_path, "fasta")
    print(f"  RC alignment written: {len(rc_records):,} sequences → {out_path}")


# ---------------------------------------------------------------------------
# 8–9. DEGEPRIME wrappers
# ---------------------------------------------------------------------------

def run_degeprime(
    degeprime_path: str,
    in_path: str,
    out_path: str,
    length: int,
    degeneracy: int,
) -> None:
    """Run DegePrime.pl for a single primer length."""
    cmd = [
        "perl", degeprime_path,
        "-i", in_path,
        "-d", str(degeneracy),
        "-l", str(length),
        "-o", out_path,
    ]
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"DegePrime.pl failed (exit {result.returncode}) for length {length}:\n"
            f"{result.stderr[-2000:]}"
        )


def run_degeprime_all_lengths(
    degeprime_path: str,
    in_path: str,
    out_dir: str,
    prefix: str,
    min_len: int,
    max_len: int,
    degeneracy: int,
) -> str:
    """
    Run DEGEPRIME for each primer length in [min_len, max_len] and merge
    all per-length TSVs into a single merged TSV (with added PrimerLen column).
    Returns path to the merged TSV.
    """
    lengths = list(range(min_len, max_len + 1))
    per_len_paths = []

    for L in lengths:
        out_path = os.path.join(out_dir, f"{prefix}_deg_{L}.tsv")
        print(f"    DEGEPRIME length {L} → {os.path.basename(out_path)}")
        run_degeprime(degeprime_path, in_path, out_path, L, degeneracy)
        per_len_paths.append((L, out_path))

    # Merge
    merged_path = os.path.join(out_dir, f"{prefix}_primers.tsv")
    with open(merged_path, "w", newline='') as out_fh:
        writer = None
        for L, path in per_len_paths:
            if not os.path.exists(path):
                continue
            with open(path, newline='') as in_fh:
                reader = csv.DictReader(in_fh, delimiter='\t')
                if writer is None:
                    fieldnames = list(reader.fieldnames or []) + ["PrimerLen"]
                    writer = csv.DictWriter(
                        out_fh, fieldnames=fieldnames, delimiter='\t',
                        extrasaction='ignore',
                    )
                    writer.writeheader()
                for row in reader:
                    row["PrimerLen"] = L
                    writer.writerow(row)

    print(f"  Merged primers → {merged_path}")
    return merged_path


# ---------------------------------------------------------------------------
# 10. Parse DEGEPRIME output
# ---------------------------------------------------------------------------

def parse_degeprime_output(tsv_path: str, min_coverage: float) -> list:
    """
    Parse a DEGEPRIME merged TSV (output of run_degeprime_all_lengths).
    Filters rows by FractionMatching >= min_coverage, NumberSpanning > 0,
    and non-trivial PrimerSeq.
    Returns list of dicts with keys: pos, seq, length, degeneracy, coverage.
    """
    candidates = []
    if not os.path.exists(tsv_path):
        return candidates

    with open(tsv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                fraction = float(row.get("FractionMatching", 0))
                spanning = int(row.get("NumberSpanning", 0))
                seq = row.get("PrimerSeq", "").strip().upper()
                pos = int(row.get("Position", 0))
                degen = int(row.get("Degeneracy", 0))
                plen = int(row.get("PrimerLen", len(seq)))
            except (ValueError, TypeError):
                continue

            if fraction < min_coverage:
                continue
            if spanning <= 0:
                continue
            if not seq or all(b in ('N', '-', '') for b in seq):
                continue

            candidates.append({
                "pos": pos,
                "seq": seq,
                "length": plen,
                "degeneracy": degen,
                "coverage": fraction,
            })

    return candidates


# ---------------------------------------------------------------------------
# 11. Convert reverse-primer positions to forward-strand coordinates
# ---------------------------------------------------------------------------

def process_reverse_primers(rev_list: list, alignment_len: int) -> list:
    """
    Convert RC-alignment positions to forward-strand coordinates.
    fwd_pos = alignment_len - rc_pos - primer_len
    The PrimerSeq from DEGEPRIME on the RC alignment is already the correct
    5'→3' reverse primer sequence (no additional RC needed).
    """
    processed = []
    for p in rev_list:
        fwd_pos = alignment_len - p["pos"] - p["length"]
        entry = dict(p)
        entry["fwd_pos"] = fwd_pos
        processed.append(entry)
    return processed


# ---------------------------------------------------------------------------
# 12. Calculate Tm
# ---------------------------------------------------------------------------

def calculate_tm(seq: str) -> float:
    """
    Calculate nearest-neighbor Tm for a primer sequence using standard PCR
    conditions: [primer]=250 nM, [Na]=50 mM, Owczarzy 2004 salt correction.
    Biopython handles IUPAC ambiguity codes.
    """
    # Strip gaps before Tm calculation
    clean = seq.replace('-', '').upper()
    return Tm_NN(clean, dnac1=250, dnac2=0, Na=50, saltcorr=5)


# ---------------------------------------------------------------------------
# 13. Find primer pairs
# ---------------------------------------------------------------------------

def find_primer_pairs(
    fwd_list: list,
    rev_list: list,
    max_amplicon_len: int,
    tm_mismatch: float,
) -> list:
    """
    Find valid primer pairs using bisect for efficiency.

    A pair (f, r) is valid if:
      1. r.fwd_pos >= f.pos + f.length      (no overlap)
      2. r.fwd_pos + r.length - f.pos <= max_amplicon_len
      3. |f.tm - r.tm| <= tm_mismatch

    Rev primers are sorted by fwd_pos; for each fwd primer we use bisect
    to find the valid window.
    """
    # Sort rev by fwd_pos for bisect
    rev_sorted = sorted(rev_list, key=lambda r: r["fwd_pos"])
    rev_positions = [r["fwd_pos"] for r in rev_sorted]

    pairs = []
    for f in fwd_list:
        f_end = f["pos"] + f["length"]
        # r.fwd_pos must be >= f_end  (no overlap)
        lo = bisect.bisect_left(rev_positions, f_end)
        # r.fwd_pos + r.length <= f.pos + max_amplicon_len
        # → r.fwd_pos <= f.pos + max_amplicon_len - r.length
        # We use a conservative upper bound; exact check below
        hi_bound = f["pos"] + max_amplicon_len
        hi = bisect.bisect_right(rev_positions, hi_bound)

        for r in rev_sorted[lo:hi]:
            amplicon_len = r["fwd_pos"] + r["length"] - f["pos"]
            if amplicon_len > max_amplicon_len:
                continue
            if abs(f["tm"] - r["tm"]) > tm_mismatch:
                continue
            pairs.append({
                "fwd_seq": f["seq"],
                "rev_seq": r["seq"],
                "fwd_pos": f["pos"],
                "rev_pos": r["fwd_pos"],
                "amplicon_len": amplicon_len,
                "fwd_tm": f["tm"],
                "rev_tm": r["tm"],
                "fwd_degeneracy": f["degeneracy"],
                "rev_degeneracy": r["degeneracy"],
                "fwd_coverage": f["coverage"],
                "rev_coverage": r["coverage"],
            })
    return pairs


# ---------------------------------------------------------------------------
# 14. Quality filters (primer dimer, homodimer, hairpin)
# ---------------------------------------------------------------------------

def _bases_can_pair(b1: str, b2: str) -> bool:
    """
    Return True if IUPAC bases b1 and b2 can form a Watson-Crick pair
    (i.e., their base sets share at least one complementary pair).
    """
    comp_b1 = {_IUPAC_COMP.get(b, 'N') for b in iupac_to_bases(b1)}
    return bool(comp_b1 & iupac_to_bases(b2))


def check_primer_dimer(seq1: str, seq2: str, window: int = 5, min_pairs: int = 4) -> bool:
    """
    Return True if a 3′-end primer dimer is likely (potential problem).
    Checks antiparallel complementarity between the last *window* bases of
    seq1 and seq2: counts positions in the antiparallel alignment that can
    form Watson-Crick pairs.
    """
    end1 = seq1[-window:].upper()
    end2 = seq2[-window:].upper()
    # Antiparallel: end2 is reversed for alignment with end1
    end2_rev = end2[::-1]
    n_pairs = sum(
        _bases_can_pair(b1, b2)
        for b1, b2 in zip(end1, end2_rev)
    )
    return n_pairs >= min_pairs


def check_hairpin(seq: str, stem: int = 4) -> bool:
    """
    Simplified 3′-end hairpin check.
    Return True if the last *stem* bases of *seq* can anneal to any upstream
    window of the same size (potential problem).
    """
    seq = seq.upper()
    if len(seq) < stem * 2 + 1:
        return False
    tail = seq[-stem:]
    tail_comp = iupac_complement(tail)[::-1]  # RC of tail
    # Check if tail_comp appears in the body (exclude last stem bases)
    body = seq[:-stem]
    for i in range(len(body) - stem + 1):
        window = body[i:i + stem]
        if all(_bases_can_pair(a, b) for a, b in zip(tail_comp, window)):
            return True
    return False


def filter_pairs(pairs: list) -> list:
    """
    Remove pairs with primer dimer, homodimer, or hairpin issues.
    """
    filtered = []
    for p in pairs:
        fwd = p["fwd_seq"]
        rev = p["rev_seq"]

        # Primer dimer (fwd 3′ vs rev 3′)
        if check_primer_dimer(fwd, rev):
            continue

        # Homodimer (fwd vs itself)
        if check_primer_dimer(fwd, fwd):
            continue

        # Homodimer (rev vs itself)
        if check_primer_dimer(rev, rev):
            continue

        # Hairpin (fwd)
        if check_hairpin(fwd):
            continue

        # Hairpin (rev)
        if check_hairpin(rev):
            continue

        filtered.append(p)
    return filtered


# ---------------------------------------------------------------------------
# 15. Write candidate pairs
# ---------------------------------------------------------------------------

PAIR_COLUMNS = [
    "fwd_seq", "rev_seq", "fwd_pos", "rev_pos", "amplicon_len",
    "fwd_tm", "rev_tm", "fwd_degeneracy", "rev_degeneracy",
    "fwd_coverage", "rev_coverage",
]


def write_candidate_pairs(pairs: list, path: str) -> None:
    """
    Write candidate primer pairs to a TSV, sorted by amplicon_len ascending
    then average coverage descending.
    """
    pairs_sorted = sorted(
        pairs,
        key=lambda p: (p["amplicon_len"], -((p["fwd_coverage"] + p["rev_coverage"]) / 2)),
    )
    with open(path, "w", newline='') as fh:
        writer = csv.DictWriter(
            fh, fieldnames=PAIR_COLUMNS, delimiter='\t', extrasaction='ignore'
        )
        writer.writeheader()
        for p in pairs_sorted:
            # Round Tm values for readability
            row = dict(p)
            row["fwd_tm"] = round(p["fwd_tm"], 2)
            row["rev_tm"] = round(p["rev_tm"], 2)
            row["fwd_coverage"] = round(p["fwd_coverage"], 4)
            row["rev_coverage"] = round(p["rev_coverage"], 4)
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Helper: get alignment length
# ---------------------------------------------------------------------------

def _alignment_length(fasta_path: str) -> int:
    """Return the alignment length (sequence length) from the first record."""
    for rec in SeqIO.parse(fasta_path, "fasta"):
        return len(rec.seq)
    return 0


# ---------------------------------------------------------------------------
# 16. Per-locus orchestrator
# ---------------------------------------------------------------------------

def process_locus(locus: str, args: argparse.Namespace) -> None:
    """Orchestrate the full primer design pipeline for one locus."""
    locus_dir = os.path.join(args.outdir, locus)
    os.makedirs(locus_dir, exist_ok=True)

    fasta_path = os.path.join(args.datadir, f"{locus}.fasta")
    coverage_path = os.path.join(args.datadir, f"{locus}_tax_coverage.txt")

    if not os.path.exists(fasta_path):
        print(f"  ERROR: Input FASTA not found: {fasta_path}")
        return

    # --- Step 1: Subsample ---
    print(f"\n  [1/9] Subsampling sequences...")
    sub_path = os.path.join(locus_dir, f"{locus}_subsampled.fasta")
    sampled = subsample_records(fasta_path, coverage_path, args.seqs_per_family, args.seed)
    if len(sampled) < 2:
        print(f"  ERROR: Only {len(sampled)} sequences after subsampling; aborting {locus}.")
        return
    SeqIO.write(sampled, sub_path, "fasta")
    print(f"  Subsampled FASTA written: {sub_path}")

    # --- Step 2: Align ---
    print(f"\n  [2/9] Aligning with MAFFT...")
    aligned_path = os.path.join(locus_dir, f"{locus}_aligned.fasta")
    run_mafft(sub_path, aligned_path, args.mafft_threads)

    # --- Step 3: Trim ---
    print(f"\n  [3/9] Trimming with TrimAl (gap threshold {args.trimal_gt})...")
    trimmed_path = os.path.join(locus_dir, f"{locus}_alignment.fasta")
    run_trimal(aligned_path, trimmed_path, args.trimal_gt)
    aln_len = _alignment_length(trimmed_path)
    print(f"  Trimmed alignment length: {aln_len:,} bp")

    if aln_len == 0:
        print(f"  ERROR: Trimmed alignment is empty; aborting {locus}.")
        return

    # --- Step 4: RC alignment ---
    print(f"\n  [4/9] Creating reverse-complement alignment...")
    rc_path = os.path.join(locus_dir, f"{locus}_rc_alignment.fasta")
    reverse_complement_alignment(trimmed_path, rc_path)

    # --- Step 5: DEGEPRIME (forward) ---
    print(f"\n  [5/9] Running DEGEPRIME on forward alignment...")
    fwd_tsv = run_degeprime_all_lengths(
        args.degeprime, trimmed_path, locus_dir,
        f"{locus}_fwd",
        args.min_primer_len, args.max_primer_len, args.max_degeneracy,
    )

    # --- Step 6: DEGEPRIME (reverse) ---
    print(f"\n  [6/9] Running DEGEPRIME on RC alignment...")
    rev_tsv_path = run_degeprime_all_lengths(
        args.degeprime, rc_path, locus_dir,
        f"{locus}_rev",
        args.min_primer_len, args.max_primer_len, args.max_degeneracy,
    )

    # --- Step 7: Parse & filter ---
    print(f"\n  [7/9] Parsing and filtering DEGEPRIME output...")
    fwd_raw = parse_degeprime_output(fwd_tsv, args.min_coverage)
    rev_raw = parse_degeprime_output(rev_tsv_path, args.min_coverage)
    print(f"  Forward candidates (pre-Tm): {len(fwd_raw):,}")
    print(f"  Reverse candidates (pre-Tm): {len(rev_raw):,}")

    if not fwd_raw or not rev_raw:
        print(
            f"  WARNING: No candidates passed coverage filter for {locus}. "
            f"Try --min-coverage 0.5 or --max-degeneracy 96."
        )
        pairs_path = os.path.join(locus_dir, f"{locus}_candidate_pairs.tsv")
        write_candidate_pairs([], pairs_path)
        return

    # --- Step 8: Calculate Tm ---
    print(f"\n  [8/9] Calculating melting temperatures...")
    for p in fwd_raw:
        p["tm"] = calculate_tm(p["seq"])
    for p in rev_raw:
        p["tm"] = calculate_tm(p["seq"])

    # Convert rev positions to forward-strand coordinates
    rev_processed = process_reverse_primers(rev_raw, aln_len)

    # Write rev primers with fwd_pos added (overwrite merged TSV)
    rev_tsv_out = os.path.join(locus_dir, f"{locus}_rev_primers.tsv")
    with open(rev_tsv_out, "w", newline='') as fh:
        fieldnames = ["pos", "fwd_pos", "seq", "length", "degeneracy", "coverage", "tm"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for r in rev_processed:
            writer.writerow(r)

    # Also write fwd primers
    fwd_tsv_out = os.path.join(locus_dir, f"{locus}_fwd_primers.tsv")
    with open(fwd_tsv_out, "w", newline='') as fh:
        fieldnames = ["pos", "seq", "length", "degeneracy", "coverage", "tm"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for f in fwd_raw:
            writer.writerow(f)

    # --- Step 9: Find & filter pairs ---
    print(f"\n  [9/9] Finding primer pairs (max amplicon {args.max_amplicon_len} bp)...")
    raw_pairs = find_primer_pairs(
        fwd_raw, rev_processed, args.max_amplicon_len, args.tm_mismatch
    )
    print(f"  Raw pairs before quality filtering: {len(raw_pairs):,}")

    filtered = filter_pairs(raw_pairs)
    print(f"  Pairs after quality filtering: {len(filtered):,}")

    # --- Write output ---
    pairs_path = os.path.join(locus_dir, f"{locus}_candidate_pairs.tsv")
    write_candidate_pairs(filtered, pairs_path)
    print(f"  Candidate pairs written → {pairs_path}")

    if not filtered:
        print(
            f"\n  WARNING: Zero pairs found for {locus} under "
            f"--max-amplicon-len {args.max_amplicon_len}.\n"
            f"  Suggestions:\n"
            f"    • Try --max-amplicon-len 100 to confirm the pipeline works.\n"
            f"    • Lower --min-coverage (current: {args.min_coverage}).\n"
            f"    • Raise --max-degeneracy (current: {args.max_degeneracy})."
        )
    else:
        best = filtered[0]
        print(
            f"\n  Top pair: amplicon {best['amplicon_len']} bp | "
            f"Tm fwd {best['fwd_tm']:.1f}°C rev {best['rev_tm']:.1f}°C | "
            f"avg coverage {(best['fwd_coverage']+best['rev_coverage'])/2:.2%}"
        )


# ---------------------------------------------------------------------------
# 17. CLI and entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Design degenerate primer pairs for anuran metabarcoding loci.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--degeprime",
        required=True,
        metavar="PATH",
        help="Path to DegePrime.pl (e.g. ./DEGEPRIME/DegePrime.pl)",
    )
    parser.add_argument(
        "--loci",
        nargs="+",
        choices=ALL_LOCI,
        default=ALL_LOCI,
        metavar="LOCUS",
        help=f"Loci to process: {' '.join(ALL_LOCI)} (default: all four)",
    )
    parser.add_argument(
        "--datadir",
        default="amphibian_locus_data",
        metavar="DIR",
        help="Directory containing input FASTA and coverage files (default: amphibian_locus_data)",
    )
    parser.add_argument(
        "--outdir",
        default="primer_design",
        metavar="DIR",
        help="Root output directory (default: primer_design)",
    )
    parser.add_argument(
        "--seqs-per-family",
        type=int,
        default=5,
        metavar="N",
        dest="seqs_per_family",
        help="Max sequences per genus group for alignment (default: 5)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for subsampling (default: 42)",
    )
    parser.add_argument(
        "--min-primer-len",
        type=int,
        default=18,
        metavar="N",
        dest="min_primer_len",
        help="Minimum primer length (default: 18)",
    )
    parser.add_argument(
        "--max-primer-len",
        type=int,
        default=22,
        metavar="N",
        dest="max_primer_len",
        help="Maximum primer length (default: 22)",
    )
    parser.add_argument(
        "--max-degeneracy",
        type=int,
        default=48,
        metavar="N",
        dest="max_degeneracy",
        help="DEGEPRIME max degeneracy — must be 2^i × 3^j (default: 48)",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        metavar="F",
        dest="min_coverage",
        help="Min FractionMatching from DEGEPRIME (default: 0.8)",
    )
    parser.add_argument(
        "--max-amplicon-len",
        type=int,
        default=65,
        metavar="N",
        dest="max_amplicon_len",
        help="Max total amplicon length in bp, including primers (default: 65)",
    )
    parser.add_argument(
        "--tm-mismatch",
        type=float,
        default=5.0,
        metavar="C",
        dest="tm_mismatch",
        help="Max Tm difference within a primer pair in °C (default: 5.0)",
    )
    parser.add_argument(
        "--trimal-gt",
        type=float,
        default=0.6,
        metavar="F",
        dest="trimal_gt",
        help="TrimAl gap threshold — fraction of non-gap required (default: 0.6)",
    )
    parser.add_argument(
        "--mafft-threads",
        type=int,
        default=4,
        metavar="N",
        dest="mafft_threads",
        help="Number of threads for MAFFT (default: 4)",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Validate degeneracy
    try:
        validate_degeneracy(args.max_degeneracy)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Check dependencies
    check_dependencies(args.degeprime)

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Output directory : {os.path.abspath(args.outdir)}/")
    print(f"Input directory  : {os.path.abspath(args.datadir)}/")
    print(f"Loci             : {', '.join(args.loci)}")
    print(f"Max amplicon     : {args.max_amplicon_len} bp")
    print(f"Seqs per genus   : {args.seqs_per_family}")
    print(f"Primer lengths   : {args.min_primer_len}–{args.max_primer_len} bp")
    print(f"Max degeneracy   : {args.max_degeneracy}")
    print(f"Min coverage     : {args.min_coverage:.0%}")

    for locus in args.loci:
        print(f"\n{'='*60}")
        print(f"  Locus : {locus}")
        print(f"{'='*60}")
        try:
            process_locus(locus, args)
        except Exception as e:
            print(f"\n  ERROR processing {locus}: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{'='*60}")
    print("All done.")


if __name__ == "__main__":
    main()
