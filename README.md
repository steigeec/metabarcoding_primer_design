# metabarcoding_primer_design

Scripts for designing metabarcoding primers targeting anuran (frog & toad) DNA.

## fetch_anuran_loci.py

Downloads anuran sequences from NCBI for four metabarcoding loci and reports taxonomic family coverage.

**Loci:** 12S rRNA · 16S rRNA · cytochrome b (cytb) · COI

**Outputs** (written to `amphibian_locus_data/`):

| File | Description |
|------|-------------|
| `<LOCUS>.fasta` | All downloaded sequences |
| `<LOCUS>_tax_coverage.txt` | Sequences per anuran family |

### Install dependencies

```bash
pip install biopython
```

### Usage

```bash
# Download all four loci
python fetch_anuran_loci.py --email your@email.com

# Download specific loci
python fetch_anuran_loci.py --email your@email.com --loci 12S COI

# Use an NCBI API key for higher rate limits (10 req/s vs 3 req/s)
python fetch_anuran_loci.py --email your@email.com --api-key YOUR_KEY

# Limit sequences per locus and change output directory
python fetch_anuran_loci.py --email your@email.com --max-records 10000 --outdir my_data
```

An NCBI API key can be obtained for free at https://www.ncbi.nlm.nih.gov/account/.
