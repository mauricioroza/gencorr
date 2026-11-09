#!/usr/bin/env python3
"""
parse_sumher.py - collect LDAK SumHer (--sum-cors) results into the gencorr
standard results table.

Inputs (one pair per LDAK run, as written by the SUMHER_RG module):
  --logs  <pair>.log    LDAK's screen output; names the two summary files
                        ('--summary <trait1>.sumher.tsv', '--summary2 ...')
  --cors  <pair>.cors   'Component Value SE' table: Her1_All, Her2_All,
                        Coher_All, Cor_All, Overlap, and Scaling1/Scaling2
                        when LDAK models multiplicative inflation. Empty when
                        LDAK failed.

Output (--out prefix): <prefix>.sumher.tsv in gencorr_common.RESULT_COLUMNS.

rg, its SE and both heritabilities are LDAK's own numbers (6 decimals). LDAK
reports no p-value, so z = rg / SE and a two-sided normal p are computed here
(as BigaGWAS shows them), and the 95% CI is rg +- 1.96 SE; extra.p_method says
so. A failed or empty LDAK run gives status 'sumher_error' with LDAK's message.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import RESULT_COLUMNS, die, validate_results, write_tsv  # noqa: E402

SUMSTATS_SUFFIX = ".sumher.tsv"


def trait_ids(log_text: str, log_name: str) -> tuple[str, str]:
    ids = []
    for flag in ("--summary", "--summary2"):
        m = re.search(rf"^{flag} (\S+)$", log_text, re.MULTILINE)
        if not m or not m.group(1).endswith(SUMSTATS_SUFFIX):
            die(f"{log_name}: no '{flag} <id>{SUMSTATS_SUFFIX}' line")
        ids.append(Path(m.group(1)).name[: -len(SUMSTATS_SUFFIX)])
    return ids[0], ids[1]


def read_cors(path: Path) -> dict[str, tuple[float, float]]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    df = pd.read_csv(path, sep=r"\s+")
    se_col = "SE" if "SE" in df.columns else "SD"     # LDAK 5 wrote 'SD'
    return {r["Component"]: (float(r["Value"]), float(r[se_col])) for _, r in df.iterrows()}


def snps_used(log_text: str, sumstats_name: str) -> float:
    """Predictors of the tagging file found in one trait's summary statistics
    ('contains summary statistics for only X of the Y predictors'); all Y when
    LDAK did not warn."""
    m = re.search(rf"{re.escape(sumstats_name)} contains summary statistics for only (\d+) of (?:the )?(\d+)", log_text)
    if m:
        return float(m.group(1))
    m = re.search(r"Reading details for (\d+) predictors from", log_text)
    return float(m.group(1)) if m else np.nan


def parse_pair(log_path: Path, cors_path: Path) -> dict:
    text = log_path.read_text(errors="replace")
    t1, t2 = trait_ids(text, log_path.name)
    cors = read_cors(cors_path)
    rg, rg_se = cors.get("Cor_All", (np.nan, np.nan))
    h2_1, h2_1_se = cors.get("Her1_All", (np.nan, np.nan))
    h2_2, h2_2_se = cors.get("Her2_All", (np.nan, np.nan))
    error = next((line.strip() for line in text.splitlines() if line.startswith("Error")), None)

    if not cors:
        status = "sumher_error"
    elif np.isnan(rg) or np.isnan(rg_se) or rg_se <= 0:
        status = "rg_not_estimable"
    else:
        status = "ok"
    z = rg / rg_se if status == "ok" else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if status == "ok" else np.nan

    extra = {"coher": cors.get("Coher_All", (np.nan,))[0],
             "coher_se": cors.get("Coher_All", (np.nan, np.nan))[1],
             "overlap": cors.get("Overlap", (np.nan,))[0],
             "overlap_se": cors.get("Overlap", (np.nan, np.nan))[1]}
    if "Scaling1" in cors:
        extra.update(scaling_1=cors["Scaling1"][0], scaling_1_se=cors["Scaling1"][1],
                     scaling_2=cors["Scaling2"][0], scaling_2_se=cors["Scaling2"][1])
    extra.update(n_snps_1=snps_used(text, f"{t1}{SUMSTATS_SUFFIX}"),
                 n_snps_2=snps_used(text, f"{t2}{SUMSTATS_SUFFIX}"),
                 p_method="normal approximation of rg / SE (LDAK reports SEs only)",
                 error=error, log=log_path.name)
    extra = {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in extra.items()}

    return {
        "method": "sumher", "scope": "global", "trait1": t1, "trait2": t2,
        "locus": np.nan, "chr": np.nan, "start": np.nan, "stop": np.nan, "n_snps": np.nan,
        "rg": rg, "rg_se": rg_se,
        "rg_ci_low": rg - 1.96 * rg_se if status == "ok" else np.nan,
        "rg_ci_high": rg + 1.96 * rg_se if status == "ok" else np.nan,
        "z": z, "p": p,
        "h2_1": h2_1, "h2_1_se": h2_1_se, "h2_2": h2_2, "h2_2_se": h2_2_se,
        "status": status, "extra": json.dumps(extra, sort_keys=True),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--logs", nargs="+", required=True, help="<pair>.log files")
    p.add_argument("--cors", nargs="*", default=[], help="<pair>.cors files")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    cors_by_stem = {Path(c).name[: -len(".cors")]: Path(c) for c in args.cors}
    rows = []
    for log in sorted(args.logs, key=lambda x: Path(x).name):
        stem = Path(log).name[: -len(".log")]
        rows.append(parse_pair(Path(log), cors_by_stem.get(stem, Path("/nonexistent"))))
    df = validate_results(pd.DataFrame(rows, columns=RESULT_COLUMNS), "sumher")
    write_tsv(df, f"{args.out}.sumher.tsv")
    print(f"sumher: {len(df)} pairs, {int((df['status'] == 'ok').sum())} ok")


if __name__ == "__main__":
    main()
