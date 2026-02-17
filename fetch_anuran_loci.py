#!/usr/bin/env python3
"""
fetch_anuran_loci.py

Download anuran (Anura) sequences from NCBI for four metabarcoding loci:
  12S rRNA, 16S rRNA, cytochrome b (cytb), and cytochrome c oxidase I (COI).

For each locus, outputs to amphibian_locus_data/:
  <LOCUS>.fasta                — sequences in FASTA format
  <LOCUS>_tax_coverage.txt     — count of sequences per anuran family

Requirements:
  pip install biopython

Usage:
  python fetch_anuran_loci.py --email your@email.com
  python fetch_anuran_loci.py --email your@email.com --loci 12S COI
  python fetch_anuran_loci.py --email your@email.com --api-key YOUR_NCBI_KEY
"""

import os
import sys
import time
import argparse
from collections import defaultdict

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTPUT_DIR = "amphibian_locus_data"
MAX_RECORDS = 50_000          # cap per locus to avoid runaway downloads
BATCH_SIZE = 500              # sequences per efetch call
REQUEST_DELAY = 0.4           # seconds between NCBI calls (~2.5 req/s; limit is 3)

LOCI_QUERIES = {
    "12S": (
        "Anura[Organism] AND "
        '("12S ribosomal RNA"[Title] OR "12S rRNA"[Title] OR "12S"[Title]) '
        "NOT patent[Filter]"
    ),
    "16S": (
        "Anura[Organism] AND "
        '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR "16S"[Title]) '
        "NOT patent[Filter]"
    ),
    "cytb": (
        "Anura[Organism] AND "
        '("cytochrome b"[Title] OR "cytb"[Title] OR "cyt b"[Title]) '
        "NOT patent[Filter]"
    ),
    "COI": (
        "Anura[Organism] AND "
        '("cytochrome c oxidase subunit I"[Title] OR '
        '"cytochrome oxidase subunit I"[Title] OR '
        '"COI"[Title] OR "COX1"[Title] OR "CO1"[Title]) '
        "NOT patent[Filter]"
    ),
}

# ---------------------------------------------------------------------------
# NCBI helper functions
# ---------------------------------------------------------------------------

def ncbi_search(query: str) -> tuple:
    """Run esearch and return (WebEnv, QueryKey, total_count, id_list)."""
    handle = Entrez.esearch(
        db="nucleotide", term=query, retmax=MAX_RECORDS, usehistory="y"
    )
    record = Entrez.read(handle)
    handle.close()
    return (
        record["WebEnv"],
        record["QueryKey"],
        int(record["Count"]),
        record["IdList"],
    )


def fetch_fasta(webenv: str, query_key: str, count: int) -> list:
    """Fetch sequences in FASTA format using NCBI server-side history."""
    sequences = []
    for start in range(0, count, BATCH_SIZE):
        batch_n = min(BATCH_SIZE, count - start)
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",
            retmode="text",
            retstart=start,
            retmax=batch_n,
            webenv=webenv,
            query_key=query_key,
        )
        batch = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        sequences.extend(batch)
        print(
            f"    {len(sequences):,}/{count:,} sequences fetched...",
            end="\r",
            flush=True,
        )
        time.sleep(REQUEST_DELAY)
    print()
    return sequences


def get_taxids(id_list: list) -> dict:
    """
    Fetch the TaxId for each sequence ID (GI or accession) via epost + esummary.
    Using epost avoids HTTP 414 errors when id_list is large.
    Returns {id_str: taxid_str}.
    """
    id_taxid: dict = {}
    total = len(id_list)
    fetched = 0
    chunk_size = 10_000  # epost handles up to ~10k IDs comfortably
    for chunk_start in range(0, total, chunk_size):
        chunk = id_list[chunk_start : chunk_start + chunk_size]
        # Post IDs to NCBI server history — no IDs in the URL
        handle = Entrez.epost(db="nucleotide", id=",".join(chunk))
        post = Entrez.read(handle)
        handle.close()
        time.sleep(REQUEST_DELAY)
        webenv = post["WebEnv"]
        query_key = post["QueryKey"]
        # Fetch summaries in batches via the server-side history
        for start in range(0, len(chunk), 200):
            handle = Entrez.esummary(
                db="nucleotide",
                webenv=webenv,
                query_key=query_key,
                retstart=start,
                retmax=200,
            )
            summaries = Entrez.read(handle)
            handle.close()
            for s in summaries:
                id_taxid[s["Id"]] = str(s["TaxId"])
            fetched += min(200, len(chunk) - start)
            print(
                f"    {fetched:,}/{total:,} taxonomy IDs fetched...",
                end="\r",
                flush=True,
            )
            time.sleep(REQUEST_DELAY)
    print()
    return id_taxid


def get_families(taxid_set: set) -> dict:
    """
    Resolve each unique TaxId to its family-level taxonomy.
    Returns {taxid_str: family_name_str}.
    """
    taxid_family: dict = {}
    taxids = sorted(taxid_set)
    total = len(taxids)
    fetched = 0
    for i in range(0, total, 200):
        batch = taxids[i : i + 200]
        handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        for rec in records:
            taxid = str(rec["TaxId"])   # str() ensures key type matches callers
            family = "Unknown"
            for node in rec.get("LineageEx", []):
                if node["Rank"] == "family":
                    family = node["ScientificName"]
            taxid_family[taxid] = family
        fetched += len(batch)
        print(
            f"    {fetched:,}/{total:,} families resolved...",
            end="\r",
            flush=True,
        )
        time.sleep(REQUEST_DELAY)
    print()
    return taxid_family


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_coverage_report(
    path: str,
    locus: str,
    family_counts: dict,
    total_seqs: int,
) -> None:
    n_families = len(family_counts)
    with open(path, "w") as fh:
        fh.write(f"Taxonomic family coverage report — {locus} (Anura)\n")
        fh.write("=" * 60 + "\n")
        fh.write(f"Total sequences  : {total_seqs:,}\n")
        fh.write(f"Families found   : {n_families}\n\n")
        fh.write(f"{'Family':<45} {'Sequences':>10}\n")
        fh.write(f"{'-'*45} {'-'*10}\n")
        for family, count in sorted(family_counts.items()):
            fh.write(f"{family:<45} {count:>10,}\n")


# ---------------------------------------------------------------------------
# Per-locus pipeline
# ---------------------------------------------------------------------------

def taxonomy_from_fasta(locus: str) -> None:
    """Recompute family coverage from an already-downloaded FASTA (no re-download)."""
    fasta_path = os.path.join(OUTPUT_DIR, f"{locus}.fasta")
    if not os.path.exists(fasta_path):
        print(f"  No FASTA found at {fasta_path} — skipping.")
        return

    print(f"  Reading {fasta_path}...")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    n = len(records)
    print(f"  {n:,} sequences loaded.")

    # Extract accession IDs from the FASTA (first token of each header)
    accessions = [rec.id for rec in records]

    print("  Fetching taxonomy IDs...")
    id_taxid = get_taxids(accessions)

    taxid_counts: dict = defaultdict(int)
    for taxid in id_taxid.values():
        taxid_counts[taxid] += 1

    unique_taxids = set(id_taxid.values())
    print(f"  Resolving {len(unique_taxids):,} unique TaxIds to families...")
    taxid_family = get_families(unique_taxids)

    family_counts: dict = defaultdict(int)
    for taxid, seq_count in taxid_counts.items():
        family = taxid_family.get(taxid, "Unknown")
        family_counts[family] += seq_count

    report_path = os.path.join(OUTPUT_DIR, f"{locus}_tax_coverage.txt")
    write_coverage_report(report_path, locus, dict(family_counts), n)
    print(f"  Wrote  : coverage report → {report_path}")
    print(f"  Families represented: {len(family_counts)}")


def process_locus(locus: str, query: str) -> None:
    print(f"  Query : {query}")

    # 1. Search
    webenv, query_key, total_count, gi_list = ncbi_search(query)
    count = min(total_count, MAX_RECORDS)
    gi_list = gi_list[:count]   # esearch returns up to retmax IDs
    print(f"  Found : {total_count:,} records (downloading {count:,})")

    if count == 0:
        print("  No records found — skipping.")
        return

    # 2. Download and write FASTA
    print("  Downloading sequences...")
    sequences = fetch_fasta(webenv, query_key, count)
    fasta_path = os.path.join(OUTPUT_DIR, f"{locus}.fasta")
    n_written = SeqIO.write(sequences, fasta_path, "fasta")
    print(f"  Wrote  : {n_written:,} sequences → {fasta_path}")

    # 3. Fetch TaxIds for all GIs returned by esearch
    print("  Fetching taxonomy IDs...")
    gi_taxid = get_taxids(gi_list)

    # 4. Count sequences per TaxId
    taxid_counts: dict = defaultdict(int)
    for taxid in gi_taxid.values():
        taxid_counts[taxid] += 1

    # 5. Resolve each unique TaxId to family
    unique_taxids = set(gi_taxid.values())
    print(f"  Resolving {len(unique_taxids):,} unique TaxIds to families...")
    taxid_family = get_families(unique_taxids)

    # 6. Aggregate to family-level counts
    family_counts: dict = defaultdict(int)
    for taxid, seq_count in taxid_counts.items():
        family = taxid_family.get(taxid, "Unknown")
        family_counts[family] += seq_count

    # 7. Write coverage report
    report_path = os.path.join(OUTPUT_DIR, f"{locus}_tax_coverage.txt")
    write_coverage_report(report_path, locus, dict(family_counts), n_written)
    print(f"  Wrote  : coverage report → {report_path}")
    print(f"  Families represented: {len(family_counts)}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    global MAX_RECORDS, OUTPUT_DIR
    parser = argparse.ArgumentParser(
        description="Download anuran sequences from NCBI for metabarcoding loci.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email address for NCBI Entrez (required by NCBI; not stored elsewhere)",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key (optional; increases rate limit from 3 to 10 req/s)",
    )
    parser.add_argument(
        "--loci",
        nargs="+",
        choices=list(LOCI_QUERIES.keys()),
        default=list(LOCI_QUERIES.keys()),
        metavar="LOCUS",
        help="Loci to download: 12S 16S cytb COI (default: all four)",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=MAX_RECORDS,
        metavar="N",
        help=f"Maximum sequences to download per locus (default: {MAX_RECORDS:,})",
    )
    parser.add_argument(
        "--outdir",
        default=OUTPUT_DIR,
        help=f"Output directory (default: {OUTPUT_DIR})",
    )
    parser.add_argument(
        "--taxonomy-only",
        action="store_true",
        help=(
            "Skip sequence download; recompute family coverage from existing "
            "FASTA files in --outdir. Useful for fixing coverage reports without "
            "re-downloading."
        ),
    )
    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    MAX_RECORDS = args.max_records
    OUTPUT_DIR = args.outdir

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {os.path.abspath(OUTPUT_DIR)}/")

    for locus in args.loci:
        print(f"\n{'='*60}")
        print(f"  Locus : {locus}")
        if args.taxonomy_only:
            taxonomy_from_fasta(locus)
        else:
            process_locus(locus, LOCI_QUERIES[locus])

    print(f"\nAll done.")


if __name__ == "__main__":
    main()
