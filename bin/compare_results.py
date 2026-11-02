#!/usr/bin/env python3
"""
compare_results.py - check that two sets of genetic-correlation results agree
within stated tolerances. This is the tool for validating gencorr against
BigaGWAS, against a direct run of LDSC/LAVA, or against a previous gencorr run.

Typical uses
------------
  # gencorr vs. the CSV BigaGWAS lets you download for an LDSC job
  compare_results.py --test results/gencorr_results.tsv --method ldsc \\
      --reference biga_ldsc.csv --reference-format bigagwas_ldsc \\
      --trait-map biga_names.csv --report ldsc_vs_biga

  # two gencorr runs (e.g. before/after a change, or chunked vs. unchunked)
  compare_results.py --test new/gencorr_results.tsv \\
      --reference old/gencorr_results.tsv --method lava

  # any other table: say which column holds what
  compare_results.py --test results/gencorr_results.tsv --method lava \\
      --reference their_lava.tsv --reference-format custom \\
      --colmap trait1=phen1,trait2=phen2,locus=locus,rg=rho,p=p

How pairs are matched
---------------------
  (trait1, trait2, locus) with the two traits ORDER-INSENSITIVE, because rg is
  symmetric. Use --trait-map (CSV: reference_name,gencorr_id) when the
  reference uses different trait names.

Tolerances (override with --tol col=value,...)
----------------------------------------------
  Absolute differences, except p which is compared as |log10(p_test) -
  log10(p_ref)|. Defaults reflect output precision: LDSC prints 4 decimals,
  LAVA p-values and CIs are Monte Carlo with 1e4-1e6 draws.

Exit status is 0 when every matched value is within tolerance and both sides
contain the same pairs (unless --allow-missing), otherwise 1. A per-pair TSV
and a short text report are written with --report.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import NA_VALUES, die, read_tsv  # noqa: E402

DEFAULT_TOL = {
    "ldsc": {"rg": 2e-4, "rg_se": 2e-4, "z": 5e-3, "p": 0.01, "h2_1": 2e-4, "h2_2": 2e-4},
    "lava": {"rg": 1e-4, "rg_ci_low": 0.02, "rg_ci_high": 0.02, "p": 0.1,
             "h2_1": 1e-6, "h2_2": 1e-6},
}

BIGAGWAS_LDSC = {  # column names in BigaGWAS' LDSC CSV export (wrap_up_massive_ldsc)
    "trait1": "Trait1", "trait2": "Trait2", "rg": "rg", "rg_se": "se", "z": "z",
    "p": "p_value",
}
BIGAGWAS_NA = ["9999", "9999.0"]


def load_reference(path, fmt, colmap):
    if fmt == "gencorr":
        return read_tsv(path)
    if fmt == "bigagwas_ldsc":
        mapping = BIGAGWAS_LDSC
        sep = ","
    elif fmt == "custom":
        if not colmap:
            die("--reference-format custom needs --colmap")
        mapping = colmap
        sep = None
    else:
        die(f"unknown reference format {fmt}")
    df = pd.read_csv(path, sep=sep, engine="python", na_values=NA_VALUES + BIGAGWAS_NA,
                     keep_default_na=False, dtype=str)
    missing = [src for src in mapping.values() if src not in df.columns]
    if missing:
        die(f"reference is missing columns {missing}; columns are {list(df.columns)}")
    out = pd.DataFrame({std: df[src] for std, src in mapping.items()})
    for col in out.columns:
        if col not in ("trait1", "trait2", "locus"):
            out[col] = pd.to_numeric(out[col], errors="coerce")
    if "locus" not in out:
        out["locus"] = np.nan
    return out


def key_frame(df: pd.DataFrame, trait_map: dict[str, str]) -> pd.DataFrame:
    df = df.copy()
    df["trait1"] = df["trait1"].astype(str).map(lambda x: trait_map.get(x, x)).astype(object)
    df["trait2"] = df["trait2"].astype(str).map(lambda x: trait_map.get(x, x)).astype(object)
    t1, t2 = df["trait1"], df["trait2"]
    a = t1.where(t1 <= t2, t2)       # order-insensitive pair key
    b = t2.where(t1 <= t2, t1)
    loc = df["locus"].astype(object).where(df["locus"].notna(), "").astype(str)
    loc = loc.replace({"nan": "", "NA": "", "<NA>": ""})
    loc = loc.str.replace(r"\.0$", "", regex=True)
    df["_key"] = a + "|" + b + "|" + loc
    return df


def compare(test: pd.DataFrame, ref: pd.DataFrame, tol: dict[str, float]):
    if test["_key"].duplicated().any():
        die("test table has duplicated pairs after matching; filter with --method")
    if ref["_key"].duplicated().any():
        die("reference table has duplicated pairs after matching")
    merged = test.merge(ref, on="_key", how="outer", suffixes=("_test", "_ref"),
                        indicator=True)
    both = merged[merged["_merge"] == "both"].copy()
    summary = {"matched": len(both),
               "only_in_test": int((merged["_merge"] == "left_only").sum()),
               "only_in_reference": int((merged["_merge"] == "right_only").sum()),
               "columns": {}}
    fail_any = pd.Series(False, index=both.index)
    for col, limit in tol.items():
        if f"{col}_test" not in both or f"{col}_ref" not in both:
            continue
        x = pd.to_numeric(both[f"{col}_test"], errors="coerce")
        y = pd.to_numeric(both[f"{col}_ref"], errors="coerce")
        if col == "p":
            diff = (np.log10(x.clip(lower=1e-300)) - np.log10(y.clip(lower=1e-300))).abs()
        else:
            diff = (x - y).abs()
        na_mismatch = x.isna() != y.isna()
        fail = (diff > limit) | na_mismatch
        both[f"{col}_diff"] = diff
        fail_any |= fail.fillna(False)
        valid = x.notna() & y.notna()
        summary["columns"][col] = {
            "tolerance": limit, "compared": int(valid.sum()),
            "max_diff": float(diff[valid].max()) if valid.any() else None,
            "n_fail": int(fail.sum()), "na_mismatch": int(na_mismatch.sum()),
            "pearson_r": float(np.corrcoef(x[valid], y[valid])[0, 1]) if valid.sum() > 2 else None,
        }
    both["within_tolerance"] = ~fail_any
    summary["n_fail"] = int(fail_any.sum())
    return merged, both, summary


def format_report(summary: dict) -> str:
    lines = [f"matched pairs      : {summary['matched']}",
             f"only in test       : {summary['only_in_test']}",
             f"only in reference  : {summary['only_in_reference']}",
             f"pairs out of tol.  : {summary['n_fail']}", "",
             f"{'column':<11}{'tol':>9}{'compared':>10}{'max_diff':>12}{'n_fail':>8}{'NA_mism':>9}{'r':>9}"]
    for col, s in summary["columns"].items():
        md = "NA" if s["max_diff"] is None else f"{s['max_diff']:.3g}"
        r = "NA" if s["pearson_r"] is None else f"{s['pearson_r']:.5f}"
        lines.append(f"{col:<11}{s['tolerance']:>9.3g}{s['compared']:>10}{md:>12}"
                     f"{s['n_fail']:>8}{s['na_mismatch']:>9}{r:>9}")
    return "\n".join(lines)


def parse_kv(text: str) -> dict[str, str]:
    out = {}
    for item in filter(None, (x.strip() for x in (text or "").split(","))):
        k, _, v = item.partition("=")
        if not v:
            die(f"expected key=value, got '{item}'")
        out[k.strip()] = v.strip()
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--test", required=True, help="gencorr results TSV")
    p.add_argument("--reference", required=True)
    p.add_argument("--method", required=True, help="method to compare (ldsc, lava, ...)")
    p.add_argument("--reference-format", default="gencorr",
                   choices=["gencorr", "bigagwas_ldsc", "custom"])
    p.add_argument("--colmap", default="", help="custom format: std=source,...")
    p.add_argument("--trait-map", help="CSV with columns reference_name,gencorr_id")
    p.add_argument("--tol", default="", help="override tolerances, e.g. rg=1e-3,p=0.2")
    p.add_argument("--only-ok", action="store_true",
                   help="compare only rows with status == ok in the test table")
    p.add_argument("--allow-missing", action="store_true",
                   help="do not fail when pairs are present on only one side")
    p.add_argument("--report", help="prefix for <prefix>.pairs.tsv and <prefix>.txt")
    args = p.parse_args(argv)

    test = read_tsv(args.test)
    test = test[test["method"] == args.method]
    if args.only_ok:
        test = test[test["status"] == "ok"]
    if test.empty:
        die(f"no rows for method '{args.method}' in {args.test}")
    ref = load_reference(args.reference, args.reference_format, parse_kv(args.colmap))
    if "method" in ref.columns:
        ref = ref[ref["method"] == args.method]
        if args.only_ok and "status" in ref.columns:
            ref = ref[ref["status"] == "ok"]

    trait_map = {}
    if args.trait_map:
        tm = pd.read_csv(args.trait_map, dtype=str)
        trait_map = dict(zip(tm.iloc[:, 0], tm.iloc[:, 1]))

    tol = dict(DEFAULT_TOL.get(args.method, {"rg": 1e-3, "p": 0.1}))
    tol.update({k: float(v) for k, v in parse_kv(args.tol).items()})

    merged, both, summary = compare(key_frame(test, trait_map), key_frame(ref, trait_map), tol)
    report = format_report(summary)
    print(report)
    if args.report:
        both.drop(columns=["_merge"]).to_csv(f"{args.report}.pairs.tsv", sep="\t",
                                             index=False, na_rep="NA")
        Path(f"{args.report}.txt").write_text(report + "\n")

    missing = summary["only_in_test"] + summary["only_in_reference"]
    ok = summary["n_fail"] == 0 and summary["matched"] > 0 and (args.allow_missing or missing == 0)
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
