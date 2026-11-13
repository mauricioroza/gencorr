#!/usr/bin/env python3
"""
parse_ldsc.py - collect LDSC h2 and rg log files into gencorr tables.

Inputs (file names carry the trait IDs, as written by the LDSC modules):
  --h2-logs   <id>.h2.log           observed-scale h2, one per trait
              <id>.h2_liab.log      liability-scale h2 (binary traits with a
                                    population prevalence only)
  --rg-logs   <batch>.rg.log        ldsc --rg trait1,traitA,traitB,... ; each
                                    log may hold several pairs

Outputs (--out prefix):
  <prefix>.ldsc.tsv             standard results table (gencorr_common.RESULT_COLUMNS)
  <prefix>.ldsc_h2.tsv          one row per trait: h2 (obs/liab), intercept, ratio, ...
                                and status: ok, h2_error (LDSC error in the
                                log) or h2_not_estimable (no h2 value printed)
  <prefix>.sample_overlap.txt   sampling-correlation matrix for LAVA, built as in
                                the LAVA sample-overlap vignette: diagonal = h2
                                intercept, off-diagonal = cross-trait gcov_int,
                                then cov2cor. Pairs that were not run, and all
                                pairs of a trait without a usable h2
                                intercept, are 0 (no correction).

Every number is taken verbatim from LDSC's own log output (LDSC prints most
values with 4 decimals; that is the precision of these tables).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import RESULT_COLUMNS, die, validate_results, write_tsv  # noqa: E402

SUMSTATS_SUFFIX = ".sumstats.gz"
NUM = r"(-?[0-9.]+(?:e[-+]?\d+)?|nan|NA|inf|-inf)"


def to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return np.nan


def trait_from_file(name: str) -> str:
    name = Path(name).name
    if not name.endswith(SUMSTATS_SUFFIX):
        die(f"unexpected sumstats name '{name}' (expected <id>{SUMSTATS_SUFFIX})")
    return name[: -len(SUMSTATS_SUFFIX)]


def find(pattern: str, text: str, n: int = 1):
    m = re.search(pattern, text)
    if not m:
        return (np.nan,) * n if n > 1 else np.nan
    vals = tuple(to_float(g) for g in m.groups())
    return vals if n > 1 else vals[0]


# --------------------------------------------------------------------------- #
# h2
# --------------------------------------------------------------------------- #
def parse_h2_log(path: Path) -> dict:
    text = path.read_text()
    out = {}
    if "Total Observed scale h2" in text:
        out["h2_obs"], out["h2_obs_se"] = find(rf"Total Observed scale h2: {NUM} \({NUM}\)", text, 2)
    if "Total Liability scale h2" in text:
        out["h2_liab"], out["h2_liab_se"] = find(rf"Total Liability scale h2: {NUM} \({NUM}\)", text, 2)
    if "Intercept: constrained" in text:
        out["intercept"], out["intercept_se"] = np.nan, np.nan
        out["intercept_constrained"] = True
    else:
        out["intercept"], out["intercept_se"] = find(rf"Intercept: {NUM} \({NUM}\)", text, 2)
    out["lambda_gc"] = find(rf"Lambda GC: {NUM}", text)
    out["mean_chi2"] = find(rf"Mean Chi\^2: {NUM}", text)
    out["ratio"], out["ratio_se"] = find(rf"Ratio: {NUM} \({NUM}\)", text, 2)
    out["n_snps"] = find(r"After merging with regression SNP LD, (\d+) SNPs remain", text)
    out["samp_prev"] = find(rf"--samp-prev {NUM}", text)
    out["pop_prev"] = find(rf"--pop-prev {NUM}", text)
    if "ERROR" in text or "Traceback" in text:
        messages = re.findall(r"^\w*(?:Error|Exception): .+$", text, re.MULTILINE)
        out["error"] = messages[-1].strip() if messages else "LDSC reported an error (see the log)"
    return out


def h2_status(row: dict) -> str:
    """ok, or why a trait has no observed-scale h2."""
    if row.get("error"):
        return "h2_error"
    h2 = row.get("h2_obs", np.nan)
    return "ok" if h2 is not None and not np.isnan(h2) else "h2_not_estimable"


def collect_h2(paths: list[Path]) -> pd.DataFrame:
    rows: dict[str, dict] = {}
    for path in paths:
        name = path.name
        if name.endswith(".h2_liab.log"):
            trait, kind = name[: -len(".h2_liab.log")], "liab"
        elif name.endswith(".h2.log"):
            trait, kind = name[: -len(".h2.log")], "obs"
        else:
            die(f"unexpected h2 log name '{name}'")
        parsed = parse_h2_log(path)
        row = rows.setdefault(trait, {"trait": trait})
        if kind == "obs":
            # The liability-only fields come from the .h2_liab.log. Skip them here
            # so the result does not depend on which log is read first (the
            # observed log has no --samp-prev/--pop-prev and would overwrite
            # them with NaN).
            liab_only = ("h2_liab", "h2_liab_se", "samp_prev", "pop_prev")
            row.update({k: v for k, v in parsed.items() if k not in liab_only})
        else:
            row["h2_liab"] = parsed.get("h2_liab", np.nan)
            row["h2_liab_se"] = parsed.get("h2_liab_se", np.nan)
            row["samp_prev"] = parsed.get("samp_prev", np.nan)
            row["pop_prev"] = parsed.get("pop_prev", np.nan)
            if parsed.get("error"):
                row["liab_error"] = parsed["error"]
    # A trait without h2 is never silent: status says why, and the task log
    # carries LDSC's message.
    for trait, row in rows.items():
        row["status"] = h2_status(row)
        for key, what in (("error", "h2"), ("liab_error", "liability-scale h2")):
            if row.get(key):
                print(f"WARNING: LDSC {what} failed for '{trait}': {row[key]}", file=sys.stderr)
        if row["status"] == "h2_not_estimable":
            print(f"WARNING: no observed-scale h2 in the LDSC log of '{trait}'", file=sys.stderr)
    cols = ["trait", "n_snps", "h2_obs", "h2_obs_se", "h2_liab", "h2_liab_se",
            "samp_prev", "pop_prev", "intercept", "intercept_se", "lambda_gc",
            "mean_chi2", "ratio", "ratio_se", "status"]
    df = pd.DataFrame(list(rows.values()))
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df[cols].sort_values("trait").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# rg
# --------------------------------------------------------------------------- #
def parse_rg_log(path: Path) -> pd.DataFrame:
    """Return one row per (p1, p2) with the summary-table values and n_snps."""
    text = path.read_text()
    m = re.search(r"--rg (\S+)", text)
    if not m:
        die(f"{path}: no '--rg' option found in the log header")
    files = m.group(1).split(",")
    t1 = trait_from_file(files[0])
    others = [trait_from_file(f) for f in files[1:]]

    # per-pair SNP counts and errors, in the order LDSC processed them
    n_snps, errors = {}, {}
    blocks = re.split(r"Computing rg for phenotype (\d+)/\d+", text)
    for idx in range(1, len(blocks), 2):
        k = int(blocks[idx]) - 2           # phenotype 2 is others[0]
        body = blocks[idx + 1]
        trait = others[k]
        n_snps[trait] = find(r"(\d+) SNPs with valid alleles", body)
        err = re.search(r"ERROR computing rg[^\n]*\n(?:.*\n)*?(\w*Error: [^\n]+)", body)
        if err:
            errors[trait] = err.group(1).strip()

    summary = re.search(r"Summary of Genetic Correlation Results\n(.*?)\n\s*\n", text, re.S)
    if not summary:
        die(f"{path}: no 'Summary of Genetic Correlation Results' table (did LDSC crash?)")
    lines = summary.group(1).strip().splitlines()
    header = lines[0].split()
    table = pd.DataFrame([ln.split() for ln in lines[1:]], columns=header)
    table["trait1"] = table["p1"].map(trait_from_file)
    table["trait2"] = table["p2"].map(trait_from_file)
    if not (table["trait1"] == t1).all():
        die(f"{path}: summary p1 does not match the first --rg file")
    for col in header[2:]:
        table[col] = table[col].map(to_float)
    table["n_snps"] = table["trait2"].map(n_snps)
    table["error"] = table["trait2"].map(errors)
    table["log"] = path.name
    return table


def build_results(rg: pd.DataFrame, h2: pd.DataFrame) -> pd.DataFrame:
    h2i = h2.set_index("trait")
    rows = []
    for _, r in rg.iterrows():
        status = "ok"
        if isinstance(r.get("error"), str):
            status = "ldsc_error"
        elif np.isnan(r["rg"]):
            status = "rg_not_estimable"  # e.g. negative or zero h2 for one trait
        extra = {
            "gcov_int": r.get("gcov_int"), "gcov_int_se": r.get("gcov_int_se"),
            "pair_h2_obs_2": r.get("h2_obs"), "pair_h2_obs_2_se": r.get("h2_obs_se"),
            "pair_h2_int_2": r.get("h2_int"), "pair_h2_int_2_se": r.get("h2_int_se"),
            "error": r.get("error") if isinstance(r.get("error"), str) else None,
            "log": r["log"],
        }
        extra = {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in extra.items()}

        def h2v(trait, col):
            return h2i.at[trait, col] if trait in h2i.index else np.nan

        rows.append({
            "method": "ldsc", "scope": "global",
            "trait1": r["trait1"], "trait2": r["trait2"],
            "locus": np.nan, "chr": np.nan, "start": np.nan, "stop": np.nan,
            "n_snps": r["n_snps"], "rg": r["rg"], "rg_se": r["se"],
            "rg_ci_low": r["rg"] - 1.959964 * r["se"], "rg_ci_high": r["rg"] + 1.959964 * r["se"],
            "z": r["z"], "p": r["p"],
            "h2_1": h2v(r["trait1"], "h2_obs"), "h2_1_se": h2v(r["trait1"], "h2_obs_se"),
            "h2_2": h2v(r["trait2"], "h2_obs"), "h2_2_se": h2v(r["trait2"], "h2_obs_se"),
            "status": status, "extra": json.dumps(extra, sort_keys=True),
        })
    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    return validate_results(df, "ldsc")


def sample_overlap_matrix(rg: pd.DataFrame, h2: pd.DataFrame) -> pd.DataFrame:
    traits = sorted(h2["trait"])
    cov = pd.DataFrame(0.0, index=traits, columns=traits)
    unusable = set()
    for t, intercept in zip(h2["trait"], h2["intercept"]):
        if np.isnan(intercept) or intercept <= 0:
            print(f"WARNING: no usable h2 intercept for '{t}'; its overlap row is "
                  "set to 0 (no correction)", file=sys.stderr)
            cov.loc[t, t] = 1.0
            unusable.add(t)
        else:
            cov.loc[t, t] = intercept
    seen: dict[tuple, list] = {}
    for _, r in rg.iterrows():
        if unusable & {r["trait1"], r["trait2"]}:
            continue   # no diagonal to scale gcov_int with: leave the pair at 0
        if r["trait1"] in cov.index and r["trait2"] in cov.index and not np.isnan(r["gcov_int"]):
            key = tuple(sorted((r["trait1"], r["trait2"])))
            seen.setdefault(key, []).append(r["gcov_int"])
    for (a, b), values in seen.items():
        if a == b:
            continue
        # LDSC can give slightly different intercepts depending on which trait is
        # p1/p2; average them rather than keeping an arbitrary one.
        cov.loc[a, b] = cov.loc[b, a] = float(np.mean(values))
    d = np.sqrt(np.diag(cov.values))
    corr = cov.values / np.outer(d, d)
    np.fill_diagonal(corr, 1.0)
    return pd.DataFrame(np.round(corr, 5), index=traits, columns=traits)


def write_overlap(mat: pd.DataFrame, path: str) -> None:
    # R read.table() format: header holds only column names; each row starts
    # with its row name. This is what LAVA's process.sample.overlap() expects.
    with open(path, "w") as fh:
        fh.write(" ".join(mat.columns) + "\n")
        for name, row in mat.iterrows():
            fh.write(name + " " + " ".join(f"{v:.5f}" for v in row.values) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h2-logs", nargs="+", required=True)
    p.add_argument("--rg-logs", nargs="*", default=[])
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    h2 = collect_h2([Path(x) for x in args.h2_logs])
    rg = (pd.concat([parse_rg_log(Path(x)) for x in sorted(args.rg_logs)], ignore_index=True)
          if args.rg_logs else pd.DataFrame(columns=["trait1", "trait2", "rg", "se", "z", "p",
                                                      "gcov_int", "n_snps", "error", "log"]))
    results = build_results(rg, h2) if len(rg) else pd.DataFrame(columns=RESULT_COLUMNS)
    results = results.sort_values(["trait1", "trait2"]).reset_index(drop=True)

    write_tsv(results, f"{args.out}.ldsc.tsv")
    write_tsv(h2, f"{args.out}.ldsc_h2.tsv")
    write_overlap(sample_overlap_matrix(rg, h2) if len(h2) else pd.DataFrame(),
                  f"{args.out}.sample_overlap.txt")
    print(f"ldsc: {len(h2)} traits, {len(results)} pairs, "
          f"status {results['status'].value_counts().to_dict() if len(results) else {}}")


if __name__ == "__main__":
    main()
