#!/usr/bin/env python3
"""
Fetch report for an assembly from NCBI. Adds three columns to output table:
    - TYPE: defines if a contig is autosomal (AUTO), sex-linked (chrX, chrY, chrZ or chrW) or the mitochondrion (chrM)
    - FEMALE_PLOIDY: ploidy setting with which to call a female individual for a specific contig 
    - MALE_PLOIDY: ploidy setting with which to call a male individual for a specific contig 

Usage:
    python fetch_reference_data.py GCA_903992535.2 output/GCA_903992535.2/contigs.tsv
"""

import json, sys, urllib.request

NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/sequence_reports?page_size=1000"

PLOIDY = {
    "AUTO":     (2, 2),
    "chrX":     (2, 1),
    "chrY":     (0, 1),
    "chrZ":     (1, 2),
    "chrW":     (1, 0),
    "chrM":     (1, 1),
    "UNPLACED": (2, 2), ## unplaced contigs are treated as autosomal
}

def classify(chr_name: str) -> str:
    if not chr_name:
        return "UNPLACED"
    c = chr_name.upper()
    if c in ("X", "CHRX", "chrX"):   return "chrX"
    if c in ("Y", "CHRY", "chrY"):   return "chrY"
    if c in ("Z", "CHRZ", "chrZ"):   return "chrZ"
    if c in ("W", "CHRW", "chrW"):   return "chrW"
    if c in ("MT", "M", "CHRMT", "MITO", "MITOCHONDRION"): return "chrM"
    return "AUTO"


def fetch(accession: str) -> list[dict]:
    rows = []
    url = NCBI_API.format(acc=accession)
    while url:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read())
        for s in data.get("reports", []):
            typ = classify(s.get("chr_name", ""))
            fp, mp = PLOIDY[typ]
            s["TYPE"] = typ
            s["FEMALE_PLOIDY"] = fp
            s["MALE_PLOIDY"] = mp
            rows.append(s)
        token = data.get("next_page_token")
        url = (NCBI_API.format(acc=accession) + f"&page_token={token}") if token else None
    return rows


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: get_vgp_contig_names.py <accession> <output.tsv>")
    accession, out = sys.argv[1], sys.argv[2]

    rows = fetch(accession)
    if not rows:
        sys.exit(f"No sequences returned for {accession}")

    # Collect all keys in order: NCBI fields first, then add the three columns
    our_cols = {"TYPE", "FEMALE_PLOIDY", "MALE_PLOIDY"}
    ncbi_cols = [k for k in rows[0] if k not in our_cols]
    cols = ncbi_cols + ["TYPE", "FEMALE_PLOIDY", "MALE_PLOIDY"]

    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
