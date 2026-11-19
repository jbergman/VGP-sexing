#!/usr/bin/env python3
"""
Download GenomeArk S3 data for a VGP assembly accession:
  - README.txt
  - sample_sheet.csv
  - callable_sites/depths/*.mosdepth.summary.txt

Usage:
    python fetch_sample_data_genark.py <accession> <output_folder>
"""

import sys, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

S3_BUCKET = "https://s3.amazonaws.com/genomeark"
S3_LIST   = "https://s3.amazonaws.com/genomeark?list-type=2&prefix={prefix}"
BASE      = S3_BUCKET + "/downstream_analyses/conservation_genomics/variant_calling/{acc}/{file}"


def s3_list_keys(prefix: str) -> list[str]:
    """Return all object keys under an S3 prefix (handles pagination)."""
    keys  = []
    token = None
    ns    = "http://s3.amazonaws.com/doc/2006-03-01/"
    while True:
        url = S3_LIST.format(prefix=prefix)
        if token:
            url += f"&continuation-token={urllib.parse.quote(token, safe='')}"
        with urllib.request.urlopen(url, timeout=30) as r:
            root = ET.fromstring(r.read())
        for c in root.findall(f"{{{ns}}}Contents"):
            key = c.findtext(f"{{{ns}}}Key", "")
            if key and not key.endswith("/"):
                keys.append(key)
        if root.findtext(f"{{{ns}}}IsTruncated", "false").lower() != "true":
            break
        token = root.findtext(f"{{{ns}}}NextContinuationToken")
    return keys


def s3_download(key: str, dest: Path):
    url = f"{S3_BUCKET}/{key}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: fetch_sample_data_genark.py <accession> <output_folder>")
    accession = sys.argv[1]
    outdir    = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)

    # README (also validates the accession exists on S3)
    readme_url = BASE.format(acc=accession, file="README.txt")
    try:
        with urllib.request.urlopen(readme_url) as r:
            readme_text = r.read().decode()
    except urllib.error.HTTPError as e:
        sys.exit(f"Failed to fetch README.txt: HTTP {e.code}")
    (outdir / "README.txt").write_text(readme_text)
    print("Downloaded README.txt")

    # sample_sheet.csv via direct key
    ss_key = (f"downstream_analyses/conservation_genomics/variant_calling"
              f"/{accession}/sample_sheet.csv")
    s3_download(ss_key, outdir / "sample_sheet.csv")
    print("Downloaded sample_sheet.csv")

    # callable_sites/depths/ — only *.mosdepth.summary.txt files
    depths_prefix = (f"downstream_analyses/conservation_genomics/variant_calling"
                     f"/{accession}/callable_sites/depths/")
    print(f"\nListing {depths_prefix} ...")
    all_keys      = s3_list_keys(depths_prefix)
    summary_keys  = [k for k in all_keys if k.endswith("summary.txt")]
    print(f"Found {len(summary_keys)} summary.txt files, downloading...")
    for key in summary_keys:
        fname = Path(key).name
        s3_download(key, outdir / "callable_sites" / "depths" / fname)
        print(f"  {fname}")


if __name__ == "__main__":
    main()
