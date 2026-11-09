#!/usr/bin/env python3
"""
export_sumstats.py - write a harmonized summary statistics file in the column
layout a specific method expects.

Each method declares its layout in FORMATS. To support a new method, add one
entry (see docs/extending.md); no other code needs to change.

A format entry has:
  columns   ordered list of (output_name, harmonized_column)
  optional  harmonized columns that are written only when they are non-missing
            for at least `min_coverage` of the variants. LDSC drops every row
            with a missing INFO/FRQ, so a half-empty optional column would
            silently discard variants; skipping it instead keeps them.
  sep       field separator
  compress  gzip the output or not

--snp-list restricts the export to the SNP IDs in a list (first column of a
whitespace-separated file, header optional, plain or gzipped), e.g. LDSC's
w_hm3.snplist. Only IDs are matched, not alleles: allele alignment is done by
each method (LAVA aligns to its LD reference). This reproduces BigaGWAS, which
filters its LAVA input to a fixed rsID list.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import HARMONIZED_COLUMNS, NA_VALUES, die, open_text  # noqa: E402

HEADER_TOKENS = {"snp", "rsid", "rs_id", "id", "snpid", "marker", "variant_id"}

FORMATS = {
    # Input for ldsc munge_sumstats.py. munge recomputes Z from P and takes the
    # sign from Z (--signed-sumstats Z,0), exactly as BigaGWAS fed it P + beta.
    "ldsc": {
        "columns": [("SNP", "SNP"), ("A1", "A1"), ("A2", "A2"), ("Z", "Z"),
                    ("P", "P"), ("N", "N")],
        "optional": [("INFO", "INFO"), ("FRQ", "EAF")],
        "min_coverage": 0.99,
        "sep": "\t",
        "compress": True,
    },
    # Input for LAVA process.input(): SNP, A1, A2, N and a Z statistic ("Z").
    "lava": {
        "columns": [("SNP", "SNP"), ("A1", "A1"), ("A2", "A2"), ("Z", "Z"),
                    ("N", "N")],
        "optional": [],
        "min_coverage": 1.0,
        "sep": "\t",
        "compress": False,  # LAVA reads with data.table::fread; plain is fastest
    },
    # Input for LDAK/SumHer --summary: Predictor A1 A2 n Z (space-separated).
    # Predictor names must match the tagging file: rsIDs by default, or chr:bp
    # (GRCh37) with --predictor-positions, as in LDAK's pre-computed taggings.
    "sumher": {
        "columns": [("Predictor", "SNP"), ("A1", "A1"), ("A2", "A2"), ("n", "N"),
                    ("Z", "Z")],
        "optional": [],
        "min_coverage": 1.0,
        "sep": " ",
        "compress": False,  # LDAK reads plain text
    },
    # Input for Popcorn fit --sfile1/--sfile2: rsid a1 a2 N Z (tab-separated).
    # Only Z is given (no beta/SE), so Popcorn uses the harmonized Z like every
    # other method (it would recompute Z from beta/SE otherwise).
    "popcorn": {
        "columns": [("rsid", "SNP"), ("a1", "A1"), ("a2", "A2"), ("N", "N"), ("Z", "Z")],
        "optional": [],
        "min_coverage": 1.0,
        "sep": "\t",
        "compress": False,
    },
}


def export(df: pd.DataFrame, fmt: str) -> pd.DataFrame:
    if fmt not in FORMATS:
        die(f"unknown format '{fmt}' (known: {', '.join(FORMATS)})")
    spec = FORMATS[fmt]
    missing = [c for c in HARMONIZED_COLUMNS if c not in df.columns]
    if missing:
        die(f"input is not a harmonized file, missing {missing}")
    out = pd.DataFrame({name: df[src] for name, src in spec["columns"]})
    for name, src in spec["optional"]:
        coverage = df[src].notna().mean() if len(df) else 0.0
        if coverage >= spec["min_coverage"]:
            out[name] = df[src]
    # Required columns must be complete; harmonize_sumstats guarantees this,
    # this check protects against hand-made inputs.
    required = [name for name, _ in spec["columns"]]
    bad = out[required].isna().any(axis=1)
    if bad.any():
        out = out.loc[~bad]
        print(f"WARNING: dropped {int(bad.sum())} rows with missing required "
              f"values for format '{fmt}'", file=sys.stderr)
    return out


def read_snp_list(path: str) -> set[str]:
    """IDs from the first column; a header line is skipped; rsIDs lower-cased
    (the harmonizer lower-cases rsIDs too, so 'RS123' and 'rs123' match)."""
    ids = set()
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            fields = line.split()
            if not fields:
                continue
            first = fields[0]
            if i == 0 and first.lower() in HEADER_TOKENS:
                continue
            ids.add(first.lower() if first[:2].lower() == "rs" else first)
    if not ids:
        die(f"--snp-list {path} contains no SNP IDs")
    return ids


def restrict(df: pd.DataFrame, ids: set[str]) -> pd.DataFrame:
    keep = df["SNP"].isin(ids)
    print(f"snp-list: kept {int(keep.sum())} of {len(df)} variants "
          f"({len(ids)} IDs in list)")
    return df.loc[keep]


def to_positions(df: pd.DataFrame, path: str) -> pd.DataFrame:
    """Rename SNP IDs to 'CHR:BP' using an rsID -> position table (PLINK .bim or a
    CHR/BP/SNP table; read like an rsID map). Variants without a position, and
    rsIDs or positions that occur more than once in the table, are dropped."""
    from harmonize_sumstats import read_rsid_map  # same bin/ directory
    table = read_rsid_map(Path(path))
    table = table.drop_duplicates(subset=["SNP", "CHR", "BP"])
    table = table[~table["SNP"].duplicated(keep=False)]
    pos = table["CHR"].astype(str) + ":" + table["BP"].astype(str)
    table = table[~pos.duplicated(keep=False)]
    lookup = pd.Series((table["CHR"].astype(str) + ":" + table["BP"].astype(str)).to_numpy(),
                       index=table["SNP"].to_numpy())
    new_ids = df["SNP"].map(lookup)
    keep = new_ids.notna()
    print(f"predictor-positions: {int(keep.sum())} of {len(df)} variants have a unique "
          f"position in {Path(path).name}; the rest are dropped")
    out = df.loc[keep].copy()
    out["SNP"] = new_ids[keep]
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="*.harmonized.tsv.gz")
    p.add_argument("--format", required=True, choices=sorted(FORMATS))
    p.add_argument("--out", required=True,
                   help="output file name ('.gz' is appended for compressed formats)")
    p.add_argument("--snp-list", help="keep only SNP IDs in this list (e.g. w_hm3.snplist)")
    p.add_argument("--predictor-positions",
                   help="name variants CHR:BP using this rsID -> position table (PLINK .bim or "
                        "CHR/BP/SNP table), e.g. for LDAK's pre-computed GRCh37 tagging files")
    args = p.parse_args(argv)

    df = pd.read_csv(args.input, sep="\t", na_values=NA_VALUES, keep_default_na=False,
                     dtype={"SNP": str, "CHR": str, "A1": str, "A2": str})
    if args.snp_list:
        df = restrict(df, read_snp_list(args.snp_list))
        if df.empty:
            die("no variants left after --snp-list: check the list uses the same "
                "IDs (rsIDs) as the summary statistics")
    if args.predictor_positions:
        df = to_positions(df, args.predictor_positions)
        if df.empty:
            die("no variants left after --predictor-positions: check the table has the "
                "summary statistics' rsIDs")
    out = export(df, args.format)
    spec = FORMATS[args.format]
    path = args.out
    if spec["compress"] and not path.endswith(".gz"):
        path += ".gz"   # the format decides compression; the name always says so
    out.to_csv(path, sep=spec["sep"], index=False, na_rep="NA",
               float_format="%.8g", compression={"method": "gzip", "mtime": 0} if spec["compress"] else None)
    print(f"{args.format}: wrote {len(out)} variants to {path}, columns {list(out.columns)}")


if __name__ == "__main__":
    main()
