#!/usr/bin/env python3
"""
Fetch NCBI BioSample metadata for all sample IDs in sample_sheet.csv
and write biosample_metadata.tsv.

Usage:
    python fetch_sample_data_ncbi.py <output_folder>

<output_folder> must already contain sample_sheet.csv
(produced by fetch_sample_data_genark.py).
"""

import csv, json, sys, time, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

ESEARCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           "?db=biosample&term={id}&retmode=json")
EFETCH  = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           "?db=biosample&id={uid}&retmode=xml")


def _urlopen_with_retry(url: str, retries: int = 5, timeout: int = 30) -> bytes:
    """Fetch url, retrying on transient 5xx errors with exponential backoff."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code >= 500 and attempt < retries - 1:
                wait = 2 ** attempt      # 1 s, 2 s, 4 s, 8 s …
                print(f"    WARNING: HTTP Error {e.code}: {e.reason}"
                      f"  (retry {attempt + 1}/{retries - 1} in {wait}s)",
                      flush=True)
                time.sleep(wait)
            else:
                raise
    raise RuntimeError("unreachable")


def fetch_biosample(sample_id: str) -> dict:
    # Step 1: resolve accession → numeric UID
    ids = json.loads(
        _urlopen_with_retry(ESEARCH.format(id=sample_id))
    ).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"sample_id": sample_id}
    time.sleep(0.3)   # stay within NCBI rate limit between esearch and efetch

    # Step 2: fetch full XML by numeric UID
    xml_text = _urlopen_with_retry(EFETCH.format(uid=ids[0]))
    root = ET.fromstring(xml_text)
    bs   = root.find("BioSample")
    if bs is None:
        return {"sample_id": sample_id}

    row = {"sample_id": sample_id, "accession": bs.get("accession", "")}

    title = bs.findtext(".//Description/Title", "").strip()
    if title:
        row["title"] = title

    org = bs.find(".//Organism")
    if org is not None:
        row["organism"] = org.get("taxonomy_name", "")
        row["tax_id"]   = org.get("taxonomy_id", "")

    for id_el in bs.findall(".//Ids/Id"):
        db  = id_el.get("db", "").strip()
        val = (id_el.text or "").strip()
        if db and val and db.lower() != "biosample":
            row[f"id_{db}"] = val

    for attr in bs.findall(".//Attribute"):
        name = attr.get("harmonized_name") or attr.get("attribute_name", "")
        name = name.strip()
        if name:
            row[name] = (attr.text or "").strip()

    return row


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: fetch_sample_data_ncbi.py <output_folder>")
    outdir = Path(sys.argv[1])

    sample_sheet = outdir / "sample_sheet.csv"
    if not sample_sheet.exists():
        sys.exit(f"ERROR: {sample_sheet} not found — run fetch_sample_data_genark.py first")

    with open(sample_sheet, newline="") as f:
        reader = csv.reader(f)
        next(reader)   # skip header
        sample_ids = list(dict.fromkeys(row[0] for row in reader if row))

    print(f"Fetching BioSample metadata for {len(sample_ids)} samples...")
    rows = []
    for i, sid in enumerate(sample_ids, 1):
        print(f"  [{i}/{len(sample_ids)}] {sid}", flush=True)
        try:
            rows.append(fetch_biosample(sid))
        except Exception as e:
            print(f"    WARNING: {e}")
            rows.append({"sample_id": sid})
        if i < len(sample_ids):
            time.sleep(0.4)

    # Build column order: sample_id first, then all other keys in encounter order
    all_keys = ["sample_id"]
    seen     = {"sample_id"}
    for row in rows:
        for k in row:
            if k not in seen:
                all_keys.append(k)
                seen.add(k)

    out = outdir / "biosample_metadata.tsv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"BioSample metadata → {out}")


if __name__ == "__main__":
    main()
