#!/usr/bin/env python3
"""
Add <sex_chr>A_RATIO and GENETIC_SEX columns to biosample_metadata.tsv
based on sex-chromosome / autosome coverage ratio.

XY system (chrX present):  ratio >= threshold → FEMALE,  < threshold → MALE
ZW system (chrZ present):  ratio >= threshold → MALE,    < threshold → FEMALE

Usage:
    python assign_genetic_sex.py <sample_dir> [--xa RATIO]

<sample_dir>   Output folder containing:
                 contigs.tsv  (from fetch_reference_data.py)
                 biosample_metadata.tsv
                 callable_sites/depths/<sample_id>.mosdepth.summary.txt

--xa RATIO     Sex-chr/autosome threshold (default: 0.75)
"""

import argparse, csv, sys
from pathlib import Path


def read_mosdepth(path: Path) -> dict[str, float]:
    result = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            chrom = row.get("chrom", "")
            try:
                result[chrom] = float(row.get("mean", ""))
            except ValueError:
                pass
    return result


def detect_name_col(contigs_path: Path, depths_dir: Path) -> str:
    """Find which column in contigs.tsv holds the chrom names used in mosdepth files."""
    summaries = list(depths_dir.glob("*.mosdepth.summary.txt"))
    if not summaries:
        sys.exit("ERROR: no *.mosdepth.summary.txt files found in " + str(depths_dir))
    mosdepth_chroms = set(read_mosdepth(summaries[0]).keys())

    with open(contigs_path, newline="") as f:
        first = next(csv.DictReader(f, delimiter="\t"), None)
    if first is None:
        sys.exit("ERROR: contigs.tsv is empty")

    for col in ("genbank_accession", "refseq_accession", "chr_name"):
        if col in first and first[col] in mosdepth_chroms:
            return col

    sys.exit(
        f"ERROR: cannot match contigs.tsv to mosdepth chrom names.\n"
        f"  mosdepth chroms (first 5): {sorted(mosdepth_chroms)[:5]}\n"
        f"  contigs.tsv columns: {list(first.keys())}"
    )


def get_reference_chroms(contigs_path: Path, name_col: str):
    """Return (sex_chr_name, sex_chr_type, longest_auto_name, longest_auto_bp)."""
    sex_chrs     = []   # (is_assembled, length, name, type)
    autos        = []
    with open(contigs_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            typ  = row.get("TYPE", "")
            name = row.get(name_col, "")
            if not name:
                continue
            is_assembled = int(row.get("role", "") == "assembled-molecule")
            try:
                bp = int(row.get("length", 0) or 0)
            except ValueError:
                bp = 0
            if typ in ("chrX", "chrZ"):
                sex_chrs.append((is_assembled, bp, name, typ))
            elif typ == "AUTO":
                try:
                    autos.append((int(row.get("length", 0) or 0), name))
                except ValueError:
                    pass
    if not sex_chrs:
        sys.exit("ERROR: no chrX or chrZ found in contigs.tsv")
    if not autos:
        sys.exit("ERROR: no AUTO chromosomes found in contigs.tsv")
    longest_bp, longest_name = max(autos)
    # prefer assembled-molecule; among ties take longest
    _, _, sex_chr, sex_chr_type = max(sex_chrs)
    return sex_chr, sex_chr_type, longest_name, longest_bp


def find_contigs_tsv(sample_dir: Path) -> Path:
    p = sample_dir / "contigs.tsv"
    if not p.exists():
        sys.exit("ERROR: contigs.tsv not found in " + str(sample_dir))
    return p


def main():
    parser = argparse.ArgumentParser(
        description="Add <sex_chr>A_RATIO and GENETIC_SEX columns to biosample_metadata.tsv."
    )
    parser.add_argument("sample_dir", type=Path,
                        help="Output folder containing contigs.tsv, biosample_metadata.tsv, "
                             "and callable_sites/depths/*.mosdepth.summary.txt")
    parser.add_argument("--xa", type=float, default=0.75, metavar="RATIO",
                        help="Sex-chr/autosome ratio threshold (default: 0.75)")
    args      = parser.parse_args()
    sample_dir = args.sample_dir
    threshold  = args.xa

    contigs_path = find_contigs_tsv(sample_dir)
    depths_dir   = sample_dir / "callable_sites" / "depths"
    metadata_tsv = sample_dir / "biosample_metadata.tsv"

    name_col                                    = detect_name_col(contigs_path, depths_dir)
    sex_chr_name, sex_chr_type, auto_name, auto_bp = get_reference_chroms(contigs_path, name_col)

    # XY: high ratio → FEMALE;  ZW: high ratio → MALE
    is_zw       = sex_chr_type == "chrZ"
    ratio_col   = "ZA_RATIO" if is_zw else "XA_RATIO"
    high_label  = "MALE"     if is_zw else "FEMALE"
    low_label   = "FEMALE"   if is_zw else "MALE"

    print(f"Contigs file  : {contigs_path.name}")
    print(f"Name column   : {name_col}")
    print(f"System        : {'ZW' if is_zw else 'XY'}")
    print(f"{sex_chr_type:<14}: {sex_chr_name}")
    print(f"Longest AUTO  : {auto_name}  ({auto_bp:,} bp)")
    print(f"Threshold     : {threshold}  (>= → {high_label}, < → {low_label})\n")

    with open(metadata_tsv, newline="") as f:
        reader     = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows       = list(reader)

    for col in (ratio_col, "GENETIC_SEX"):
        if col not in fieldnames:
            fieldnames.append(col)

    for row in rows:
        sid     = row.get("sample_id", "")
        summary = depths_dir / f"{sid}.mosdepth.summary.txt"
        if not summary.exists():
            row[ratio_col]     = ""
            row["GENETIC_SEX"] = "MISSING"
            print(f"  {sid}: summary file not found")
            continue
        covs        = read_mosdepth(summary)
        sex_chr_cov = covs.get(sex_chr_name)
        auto_cov    = covs.get(auto_name)
        if sex_chr_cov is None or auto_cov is None:
            row[ratio_col]     = ""
            row["GENETIC_SEX"] = "UNKNOWN"
            missing = sex_chr_name if sex_chr_cov is None else auto_name
            print(f"  {sid}: {missing} not found in summary → UNKNOWN")
            continue
        ratio              = sex_chr_cov / auto_cov if auto_cov else 0.0
        row[ratio_col]     = f"{ratio:.4f}"
        row["GENETIC_SEX"] = high_label if ratio >= threshold else low_label
        print(f"  {sid}: {sex_chr_type}={sex_chr_cov:.2f} / AUTO={auto_cov:.2f} = {ratio:.3f} → {row['GENETIC_SEX']}")

    with open(metadata_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nUpdated {metadata_tsv}")


if __name__ == "__main__":
    main()
