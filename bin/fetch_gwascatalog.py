#!/usr/bin/env python3
"""
Download a study's harmonised summary statistics from the GWAS Catalog.

Given an accession (e.g. GCST90000016), finds the study's harmonised file on
the GWAS Catalog FTP, downloads it with its metadata (<file>-meta.yaml),
checks the md5 given in the metadata, and writes a JSON file that
harmonize_sumstats.py reads with --source-info:

  accession, url, file, md5, genome_assembly, metadata date, trait, the
  samples (ancestry categories, size, case/control), n (total sample size,
  only when the metadata lists one sample group, of the --ancestry), and
  warnings.

Finding the file. Studies live in <base>/<range>/<accession>/harmonised/,
where <range> is the block of 1,000 accessions, e.g.
GCST90000001-GCST90001000/GCST90000016. Some studies have two harmonised
files (an older '<pmid>-<accession>-<efo>.h.tsv.gz' and a newer
'<accession>.h.tsv.gz'); the newer one is used. If the folder is not where
expected, the Catalog's harmonised_list.txt is searched for the exact
accession (a substring match could pick another study).

Only the standard library is used (the metadata is read with a small reader
for the Catalog's fixed YAML layout), so this runs in the gencorr-tools image.

--base can be a local folder laid out like the FTP (used by the tests).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

DEFAULT_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
ACCESSION = re.compile(r"^GCST(\d+)$")
# samplesheet ancestry -> GWAS Catalog ancestry category
CATEGORY = {"EUR": "European", "EAS": "East Asian"}
USER_AGENT = "gencorr-pipeline (https://github.com/mauricioroza/gencorr)"
RETRIES = 3


def die(msg: str) -> None:
    sys.exit(f"fetch_gwascatalog.py: ERROR: {msg}")


# --------------------------------------------------------------------------
# locating the study
# --------------------------------------------------------------------------
def study_dir(accession: str) -> str:
    """'<range>/<accession>', e.g. GCST90000001-GCST90001000/GCST90000016."""
    m = ACCESSION.match(accession)
    if not m:
        die(f"'{accession}' is not a GWAS Catalog accession (GCST followed by digits)")
    digits = m.group(1)
    n, width = int(digits), len(digits)
    lo = (n - 1) // 1000 * 1000 + 1
    hi = lo + 999
    return f"GCST{lo:0{width}d}-GCST{hi:0{width}d}/{accession}"


def is_local(base: str) -> bool:
    return not re.match(r"^[a-z]+://", base) or base.startswith("file://")


def local_path(base: str, rel: str = "") -> Path:
    root = base[len("file://"):] if base.startswith("file://") else base
    return Path(root) / rel


def open_url(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(req, timeout=timeout)


def read_text(base: str, rel: str) -> str | None:
    """Text of <base>/<rel>, or None if it does not exist."""
    if is_local(base):
        p = local_path(base, rel)
        return p.read_text() if p.is_file() else None
    for attempt in range(RETRIES):
        try:
            with open_url(f"{base}/{rel}") as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == RETRIES - 1:
                raise
        except urllib.error.URLError:
            if attempt == RETRIES - 1:
                raise
        time.sleep(5 * (attempt + 1))
    return None


def list_harmonised(base: str, rel_dir: str) -> list[str] | None:
    """Harmonised data files in <base>/<rel_dir>, or None if the folder is absent."""
    if is_local(base):
        p = local_path(base, rel_dir)
        return sorted(f.name for f in p.iterdir()) if p.is_dir() else None
    html = read_text(base, rel_dir + "/")
    if html is None:
        return None
    return sorted(set(re.findall(r'href="([^"/?]+)"', html)))


def pick_file(accession: str, names: list[str]) -> str | None:
    data = [n for n in names if n.endswith(".h.tsv.gz")]
    if f"{accession}.h.tsv.gz" in data:
        return f"{accession}.h.tsv.gz"
    if len(data) == 1:
        return data[0]
    if len(data) > 1:
        die(f"{accession}: several harmonised files ({', '.join(data)}); "
            "give the URL of the one you want in the samplesheet instead")
    return None


def locate(accession: str, base: str) -> str:
    """Path of the harmonised file relative to base."""
    rel_dir = f"{study_dir(accession)}/harmonised"
    names = list_harmonised(base, rel_dir)
    name = pick_file(accession, names) if names else None
    if name:
        return f"{rel_dir}/{name}"
    # fall back to the Catalog's own list, matching the accession exactly
    listing = read_text(base, "harmonised_list.txt") or ""
    hits = [line.strip()[2:] for line in listing.splitlines()
            if line.strip().split("/")[2:3] == [accession] and line.strip().endswith(".h.tsv.gz")]
    by_name = {h.split("/")[-1]: h for h in hits}
    name = pick_file(accession, list(by_name))
    if name:
        return by_name[name]
    die(f"{accession}: no harmonised summary statistics in the GWAS Catalog "
        f"(looked in {base}/{rel_dir}/ and harmonised_list.txt). Only harmonised "
        "studies can be given by accession; otherwise download the file and give its path.")
    return ""


# --------------------------------------------------------------------------
# metadata (<file>-meta.yaml)
# --------------------------------------------------------------------------
def parse_meta_yaml(text: str) -> dict:
    """Read the fields gencorr uses from a GWAS Catalog metadata file.

    The files are machine-written with a fixed layout: top-level 'key: value'
    lines, lists of '- item' lines, and a 'samples:' list of mappings whose
    entries start with '- key: value' (indented or not). Long values may
    continue on following lines; those lines are skipped.
    """
    meta: dict = {"samples": []}
    key = None               # current top-level key
    entry = None             # current samples entry
    sub = None               # current list key inside a samples entry
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        body = line.strip()
        top = re.match(r"^([A-Za-z_][\w-]*):\s*(.*)$", line)
        if indent == 0 and top:
            key, value = top.group(1), top.group(2).strip()
            entry = sub = None
            if value and not value.startswith(('"', "'")):
                meta[key] = value
            elif value.startswith(("'", '"')) and len(value) >= 2 and value[-1] == value[0]:
                meta[key] = value[1:-1]
            elif not value:
                meta.setdefault(key, None)
            continue
        if key == "samples":
            item = re.match(r"^-\s+(.*)$", body)
            kv = re.match(r"^([A-Za-z_][\w-]*):\s*(.*)$", item.group(1) if item else body)
            if item and kv:                       # '- key: value' starts an entry
                entry = {}
                meta["samples"].append(entry)
            if kv and entry is not None:
                sub = kv.group(1)
                val = kv.group(2).strip()
                entry[sub] = val if val else []
            elif item and entry is not None and sub is not None and isinstance(entry.get(sub), list):
                entry[sub].append(item.group(1).strip())
            continue
        item = re.match(r"^-\s+(.*)$", body)
        if item and key:
            if not isinstance(meta.get(key), list):
                meta[key] = []
            meta[key].append(item.group(1).strip())
        # anything else continues a long value: skipped
    samples = []
    for e in meta["samples"]:
        size = e.get("sample_size")
        cats = e.get("sample_ancestry_category", [])
        samples.append({
            "ancestry": cats if isinstance(cats, list) else [cats],
            "size": int(float(size)) if size not in (None, "", []) and re.match(r"^[\d.]+$", str(size)) else None,
            "case_control": str(e.get("case_control_study", "")).lower() == "true",
        })
    meta["samples"] = samples
    return meta


def total_n(samples: list[dict], ancestry: str) -> tuple[int | None, str | None, list[str]]:
    """(n, reason when n is None, warnings) for the samplesheet ancestry."""
    warnings = []
    want = CATEGORY.get(ancestry, ancestry)
    cats = sorted({c for s in samples for c in s["ancestry"]})
    if samples and want not in cats:
        warnings.append(f"the study's ancestry ({', '.join(cats) or 'not given'}) does not include "
                        f"{want}, but the samplesheet says {ancestry}")
    if not samples or any(s["size"] is None for s in samples):
        return None, "the GWAS Catalog metadata gives no sample size", warnings
    parts = "; ".join(f"{'/'.join(s['ancestry'])} {s['size']}" for s in samples)
    if cats != [want]:
        return None, (f"multi-ancestry or other-ancestry study ({parts}): the total N would not "
                      "match the LD references, so set n in the samplesheet"), warnings
    if len(samples) > 1:
        # several entries can be stages (discovery, replication) or the same sample
        # listed twice; their sum is not reliably the N of the summary statistics
        return None, (f"the study lists {len(samples)} sample groups ({parts}), which may be "
                      "stages or repeated entries rather than parts of one analysis, so set n "
                      "in the samplesheet"), warnings
    return samples[0]["size"], None, warnings


# --------------------------------------------------------------------------
# download
# --------------------------------------------------------------------------
def download(base: str, rel: str, dest: Path) -> str:
    """Copy <base>/<rel> to dest; return its md5."""
    for attempt in range(RETRIES):
        h = hashlib.md5()
        try:
            if is_local(base):
                src = local_path(base, rel)
                if not src.is_file():
                    die(f"{src} not found")
                with open(src, "rb") as fi, open(dest, "wb") as fo:
                    for block in iter(lambda: fi.read(1 << 20), b""):
                        h.update(block)
                        fo.write(block)
            else:
                with open_url(f"{base}/{rel}", timeout=300) as r, open(dest, "wb") as fo:
                    for block in iter(lambda: r.read(1 << 20), b""):
                        h.update(block)
                        fo.write(block)
            return h.hexdigest()
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            if attempt == RETRIES - 1:
                die(f"download of {base}/{rel} failed: {e}")
            print(f"download failed ({e}); retrying", file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    return ""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--accession", required=True, help="GWAS Catalog accession, e.g. GCST90000016")
    p.add_argument("--id", required=True, help="trait ID (output names)")
    p.add_argument("--ancestry", default="EUR", help="samplesheet ancestry (EUR, EAS)")
    p.add_argument("--base", default=DEFAULT_BASE, help="GWAS Catalog summary statistics root (URL or folder)")
    args = p.parse_args(argv)
    base = args.base.rstrip("/")

    rel = locate(args.accession, base)
    name = rel.split("/")[-1]
    meta_text = read_text(base, rel + "-meta.yaml")
    meta = parse_meta_yaml(meta_text) if meta_text else {"samples": []}

    md5 = download(base, rel, Path(name))
    expected = meta.get("data_file_md5sum")
    if expected and md5 != expected:
        die(f"{name}: md5 {md5} differs from the GWAS Catalog's {expected} (incomplete download?)")

    n, n_reason, warnings = total_n(meta["samples"], args.ancestry)
    if meta_text is None:
        warnings.append("no metadata file (-meta.yaml): sample size and md5 not available")
    info = {
        "source": "gwas_catalog",
        "accession": args.accession,
        "url": f"{base}/{rel}" if not is_local(base) else str(local_path(base, rel)),
        "file": name,
        "md5": md5,
        "md5_checked": bool(expected),
        "genome_assembly": meta.get("genome_assembly"),
        "file_type": meta.get("file_type"),
        "date_metadata_last_modified": meta.get("date_metadata_last_modified"),
        "trait": "; ".join(meta.get("trait_description") or []) if isinstance(meta.get("trait_description"), list)
                 else meta.get("trait_description"),
        "samples": meta["samples"],
        "case_control": any(s["case_control"] for s in meta["samples"]),
        "n": n,
        "n_reason": n_reason,
        "warnings": warnings,
    }
    Path(f"{args.id}.gwascatalog.json").write_text(json.dumps(info, indent=2) + "\n")
    for w in warnings:
        print(f"WARNING {args.id} ({args.accession}): {w}", file=sys.stderr)
    print(f"{args.id}: {args.accession} -> {name} (md5 {'checked' if expected else 'not checked'}); "
          f"n = {n if n is not None else 'not set (' + n_reason + ')'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
