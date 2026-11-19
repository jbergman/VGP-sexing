# vgp-sex-chromosome

A Snakemake workflow for assigning genetic sex to samples from the
[Vertebrate Genomes Project](https://vertebrategenomesproject.org/) (VGP)
Conservation Genomics dataset.

For each assembly accession the workflow:

1. downloads the chromosome-level sequence report from NCBI Datasets to find
   the sex chromosome and its length;
2. downloads the sample sheet, README, and per-sample mosdepth coverage
   summaries from [GenomeArk](https://genomeark.github.io/);
3. fetches BioSample metadata (sex, organism, collection site) from NCBI
   eutils; and
4. computes a sex-chromosome / autosome coverage ratio and assigns
   `MALE`, `FEMALE`, or `UNKNOWN` to every sample.

Both XY (mammals) and ZW (birds, reptiles) sex-determination systems are
supported.

## Pipeline

```
rule all
  └── assign_genetic_sex      (contigs.tsv + biosample_metadata.tsv → .sex_assigned)
        ├── fetch_sample_data_ncbi   (sample_sheet.csv → biosample_metadata.tsv)
        │     └── fetch_sample_data_genark  (contigs.tsv → sample_sheet.csv)
        └── fetch_reference_data     (accession → contigs.tsv)
```

| Script | What it does |
|---|---|
| `scripts/fetch_reference_data.py` | Downloads NCBI sequence report; writes `contigs.tsv` with chromosome lengths and roles |
| `scripts/fetch_sample_data_genark.py` | Lists the GenomeArk S3 folder; downloads `sample_sheet.csv`, `README.txt`, and mosdepth summary files |
| `scripts/fetch_sample_data_ncbi.py` | Calls NCBI eutils esearch + efetch for each BioSample ID in `sample_sheet.csv`; writes `biosample_metadata.tsv` |
| `scripts/assign_genetic_sex.py` | Reads per-sample mosdepth coverage, computes XA/ZA ratio, writes `GENETIC_SEX` column into `biosample_metadata.tsv` |

## Usage

```bash
# Dry-run (show what Snakemake would do, without executing):
snakemake -n

# Run all accessions with 4 parallel jobs:
snakemake --cores 4

# Re-run a specific accession (e.g. after a failed NCBI fetch):
snakemake --cores 1 --forcerun fetch_sample_data_ncbi \
    --config accessions="GCA_040182565.1"

# Run only specific accessions:
snakemake --cores 4 --config accessions="GCA_903992535.2,GCA_014839835.1"

# Pass accessions from a text file (one per line):
snakemake --cores 4 --config accessions="my_accessions.txt"

# Adjust the sex-chromosome / autosome ratio threshold (default 0.75):
snakemake --cores 4 --config xa=0.6
```

### Output

Each accession gets its own output folder:

```
output/
  GCA_903992535.2/
    contigs.tsv                  # chromosome lengths from NCBI
    sample_sheet.csv             # GenomeArk sample list
    biosample_metadata.tsv       # NCBI metadata + GENETIC_SEX column
    .sex_assigned                # sentinel: pipeline finished
    logs/
      fetch_reference_data.log
      fetch_sample_data_genark.log
      fetch_sample_data_ncbi.log
      assign_genetic_sex.log
```

`biosample_metadata.tsv` is the primary result. It contains one row per
sample with columns including `BIOSAMPLE_ID`, `ORGANISM`, `COLLECTED_BY`,
`SEX` (as reported to NCBI), `XA_RATIO` or `ZA_RATIO`, and `GENETIC_SEX`
(pipeline-assigned).

## Configuration

Options can be set in `config.yaml` or passed directly with `--config`:

| Key | Default | Description |
|---|---|---|
| `output_dir` | `output` | Root folder for per-accession output |
| `xa` | `0.75` | Sex-chr / autosome coverage ratio threshold |
| `accessions` | *(all from S3)* | Comma-separated list or path to text file |

## Install

```bash
conda env create -f environment.yaml
conda activate vgp_sex_chr
```

The workflow uses only Python standard library plus `pysam` (for BAM QC).
All network calls to NCBI and S3 use the standard `urllib` module — no boto3,
no biopython.

NCBI eutils calls include automatic retry with exponential back-off on HTTP
5xx errors (up to 5 attempts: 1 s, 2 s, 4 s, 8 s, 16 s).

## Notes on sex-chromosome detection

For assemblies where multiple sequence report entries match chrX or chrZ
(e.g. the main chromosome plus unlocalized scaffolds), the script always
prefers the `assembled-molecule` entry and falls back to other roles only
when no assembled molecule is annotated.

Assemblies with no chrX or chrZ in the NCBI sequence report (e.g. species
where sex chromosomes are not yet annotated) are logged as errors and
skipped — this is expected for a subset of VGP accessions.
