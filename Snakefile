"""
VGP Sex-Chromosome Workflow — Snakemake version
================================================

Usage:
    # Dry-run (show what would be done):
    snakemake -n

    # Run with 4 parallel jobs:
    snakemake --cores 4

    # Re-run everything even if outputs exist:
    snakemake --cores 4 --forceall

    # Custom sex-chr ratio threshold:
    snakemake --cores 4 --config xa=0.6

    # Limit to specific accessions:
    snakemake --cores 4 --config accessions="GCA_903992535.2,GCA_014839835.1"

Configuration (config.yaml or --config key=value):
    output_dir   Root folder for per-accession output   (default: output)
    xa           Sex-chr/autosome ratio threshold        (default: 0.75)
    accessions   Comma-separated list of accessions, or
                 path to a text file (one per line).     (default: list from S3)
"""

import urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

# ── configuration ─────────────────────────────────────────────────────────────
OUTPUT_DIR = Path(config.get("output_dir", "output"))
XA         = config.get("xa", 0.75)
SCRIPTS    = Path(workflow.basedir) / "scripts"

S3_LIST = ("https://s3.amazonaws.com/genomeark"
           "?list-type=2&prefix={prefix}&delimiter=/")
S3_PREFIX = "downstream_analyses/conservation_genomics/variant_calling/"


def list_s3_accessions():
    accessions, token = [], None
    ns = "http://s3.amazonaws.com/doc/2006-03-01/"
    while True:
        url = S3_LIST.format(prefix=S3_PREFIX)
        if token:
            url += f"&continuation-token={urllib.parse.quote(token, safe='')}"
        with urllib.request.urlopen(url, timeout=30) as r:
            root = ET.fromstring(r.read())
        for cp in root.findall(f"{{{ns}}}CommonPrefixes"):
            folder = cp.findtext(f"{{{ns}}}Prefix", "").rstrip("/").split("/")[-1]
            if folder.startswith("GCA") or folder.startswith("GCF"):
                accessions.append(folder)
        if root.findtext(f"{{{ns}}}IsTruncated", "false").lower() != "true":
            break
        token = root.findtext(f"{{{ns}}}NextContinuationToken")
    return sorted(accessions)


# resolve the accession list once at parse time
_acc_config = config.get("accessions", None)
if _acc_config:
    _raw = str(_acc_config)
    if Path(_raw).is_file():                        # path to a text file
        ACCESSIONS = [l.strip() for l in Path(_raw).read_text().splitlines() if l.strip()]
    else:                                           # comma-separated string
        ACCESSIONS = [a.strip() for a in _raw.split(",") if a.strip()]
else:
    print("Fetching accession list from GenomeArk S3 …", flush=True)
    ACCESSIONS = list_s3_accessions()
    print(f"Found {len(ACCESSIONS)} accession(s).", flush=True)


# ── target rule ───────────────────────────────────────────────────────────────
rule all:
    """Collect the final output of every accession."""
    input:
        expand(str(OUTPUT_DIR / "{acc}" / ".sex_assigned"), acc=ACCESSIONS)


# ── rule 1: fetch reference contig data from NCBI ─────────────────────────────
rule fetch_reference_data:
    """Download sequence reports from NCBI and write contigs.tsv."""
    output:
        contigs = str(OUTPUT_DIR / "{acc}" / "contigs.tsv")
    log:
        str(OUTPUT_DIR / "{acc}" / "logs" / "fetch_reference_data.log")
    params:
        script = str(SCRIPTS / "fetch_reference_data.py")
    shell:
        """
        mkdir -p $(dirname {log})
        python {params.script} {wildcards.acc} {output.contigs} \
            > {log} 2>&1
        """


# ── rule 2a: fetch GenomeArk S3 data ─────────────────────────────────────────
rule fetch_sample_data_genark:
    """Download sample_sheet.csv, README, and mosdepth summary files from S3."""
    input:
        contigs = str(OUTPUT_DIR / "{acc}" / "contigs.tsv")   # wait for rule 1
    output:
        sample_sheet = str(OUTPUT_DIR / "{acc}" / "sample_sheet.csv")
    log:
        str(OUTPUT_DIR / "{acc}" / "logs" / "fetch_sample_data_genark.log")
    params:
        script = str(SCRIPTS / "fetch_sample_data_genark.py"),
        outdir = str(OUTPUT_DIR / "{acc}")
    shell:
        """
        mkdir -p $(dirname {log})
        python {params.script} {wildcards.acc} {params.outdir} \
            > {log} 2>&1
        """


# ── rule 2b: fetch NCBI BioSample metadata ────────────────────────────────────
rule fetch_sample_data_ncbi:
    """Fetch BioSample metadata from NCBI for all samples in sample_sheet.csv."""
    input:
        sample_sheet = str(OUTPUT_DIR / "{acc}" / "sample_sheet.csv")
    output:
        metadata = str(OUTPUT_DIR / "{acc}" / "biosample_metadata.tsv")
    log:
        str(OUTPUT_DIR / "{acc}" / "logs" / "fetch_sample_data_ncbi.log")
    params:
        script = str(SCRIPTS / "fetch_sample_data_ncbi.py"),
        outdir = str(OUTPUT_DIR / "{acc}")
    shell:
        """
        python {params.script} {params.outdir} \
            > {log} 2>&1
        """


# ── rule 3: assign genetic sex ────────────────────────────────────────────────
rule assign_genetic_sex:
    """Add XA_RATIO / ZA_RATIO and GENETIC_SEX columns to biosample_metadata.tsv."""
    input:
        contigs  = str(OUTPUT_DIR / "{acc}" / "contigs.tsv"),
        metadata = str(OUTPUT_DIR / "{acc}" / "biosample_metadata.tsv")
    output:
        # assign_genetic_sex.py edits biosample_metadata.tsv in place;
        # we use a sentinel file to track completion
        done = str(OUTPUT_DIR / "{acc}" / ".sex_assigned")
    log:
        str(OUTPUT_DIR / "{acc}" / "logs" / "assign_genetic_sex.log")
    params:
        script   = str(SCRIPTS / "assign_genetic_sex.py"),
        sample_dir = str(OUTPUT_DIR / "{acc}"),
        xa       = XA
    shell:
        """
        python {params.script} {params.sample_dir} --xa {params.xa} \
            > {log} 2>&1 \
            && touch {output.done}
        """
