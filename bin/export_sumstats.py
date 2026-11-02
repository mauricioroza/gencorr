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
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import HARMONIZED_COLUMNS, NA_VALUES, die  # noqa: E402

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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="*.harmonized.tsv.gz")
    p.add_argument("--format", required=True, choices=sorted(FORMATS))
    p.add_argument("--out", required=True,
                   help="output file name ('.gz' is appended for compressed formats)")
    args = p.parse_args(argv)

    df = pd.read_csv(args.input, sep="\t", na_values=NA_VALUES, keep_default_na=False,
                     dtype={"SNP": str, "CHR": str, "A1": str, "A2": str})
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
