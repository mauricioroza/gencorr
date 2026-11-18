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

--adjust-p (used for the LAVA merge) adds multiple-testing corrections to
`extra` of every status-ok row with a p-value, per method, over all those rows
in the merged table (every pair x locus tested in the run, as in the LAVA
paper): `p_bonferroni` = min(1, p * n_tests), `q_fdr` = Benjamini-Hochberg
q-value, and `n_tests`. Rows that were not tested (univ_filtered, ...) are not
counted. The keys are appended to the JSON text, so existing values keep their
exact formatting.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
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


ADJUST_KEYS = ("p_bonferroni", "q_fdr", "n_tests")
_OLD_KEYS = re.compile(r'"(?:%s)"\s*:\s*[^,}]*,?' % "|".join(ADJUST_KEYS))


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (same as R's p.adjust(p, "BH"))."""
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    scaled = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def _add_keys(extra, values: dict) -> str:
    text = extra if isinstance(extra, str) and extra.strip().startswith("{") else "{}"
    text = _OLD_KEYS.sub("", text).rstrip()
    text = re.sub(r",\s*}$", "}", text)                   # a removed key may leave ',}'
    body = ",".join(f'"{k}":{v:.6g}' if isinstance(v, float) else f'"{k}":{v}'
                    for k, v in values.items())
    return text[:-1] + ("," if text.strip() != "{}" else "") + body + "}"


def adjust_p(df: pd.DataFrame) -> pd.DataFrame:
    """Add p_bonferroni / q_fdr / n_tests to `extra` of tested rows (see module doc)."""
    df = df.copy()
    tested = (df["status"] == "ok") & np.isfinite(pd.to_numeric(df["p"], errors="coerce"))
    for _method, idx in df[tested].groupby("method").groups.items():
        p = df.loc[idx, "p"].astype(float).to_numpy()
        m = len(p)
        bonf = np.minimum(p * m, 1.0)
        q = bh_qvalues(p)
        df.loc[idx, "extra"] = [
            _add_keys(e, {"p_bonferroni": float(b), "q_fdr": float(qq), "n_tests": m})
            for e, b, qq in zip(df.loc[idx, "extra"], bonf, q)]
    return df


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--adjust-p", action="store_true",
                   help="add p_bonferroni, q_fdr and n_tests to extra of status-ok rows")
    p.add_argument("inputs", nargs="*", help="results tables to merge")
    args = p.parse_args(argv)
    merged = merge([Path(x) for x in args.inputs])
    if args.adjust_p and len(merged):
        merged = adjust_p(merged)
    write_tsv(merged, args.out)
    counts = merged.groupby(["method", "status"]).size().to_dict() if len(merged) else {}
    print(f"merged {len(args.inputs)} tables -> {len(merged)} rows; (method, status): {counts}")


if __name__ == "__main__":
    main()
