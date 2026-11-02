#!/usr/bin/env python3
"""
merge_results.py - concatenate standard results tables (from any methods and
any chunks), validate them against the gencorr contract, sort, and write one
table.

Used twice in the pipeline: to merge LAVA locus chunks into one LAVA table,
and to merge all methods into results/gencorr_results.tsv.

Checks performed (any failure stops the pipeline with a clear message):
  * every input has exactly the RESULT_COLUMNS (see gencorr_common.py)
  * no duplicated (method, trait1, trait2, locus) rows. A duplicate means the
    same estimate was produced twice, i.e. a pairing or chunking bug.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import (RESULT_COLUMNS, die, empty_results, read_tsv,  # noqa: E402
                            validate_results, write_tsv)

KEY = ["method", "trait1", "trait2", "locus"]


def merge(paths: list[Path]) -> pd.DataFrame:
    tables = []
    for path in paths:
        df = read_tsv(path)
        tables.append(validate_results(df, str(path)))
    if not tables:
        return empty_results()
    merged = pd.concat(tables, ignore_index=True)
    dup = merged.duplicated(KEY, keep=False)
    if dup.any():
        die("duplicated estimates for (method, trait1, trait2, locus):\n"
            + merged.loc[dup, KEY].head(20).to_string(index=False))
    # natural ordering: method, traits, then chromosome/position for local results
    merged["_chr"] = pd.to_numeric(merged["chr"].astype(str).str.replace("X", "23"),
                                   errors="coerce")
    merged = merged.sort_values(["method", "trait1", "trait2", "_chr", "start"],
                                na_position="first", kind="mergesort")
    return merged.drop(columns="_chr")[RESULT_COLUMNS].reset_index(drop=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("inputs", nargs="*", help="results tables to merge")
    args = p.parse_args(argv)
    merged = merge([Path(x) for x in args.inputs])
    write_tsv(merged, args.out)
    counts = merged.groupby(["method", "status"]).size().to_dict() if len(merged) else {}
    print(f"merged {len(args.inputs)} tables -> {len(merged)} rows; (method, status): {counts}")


if __name__ == "__main__":
    main()
