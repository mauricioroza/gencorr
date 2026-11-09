#!/usr/bin/env python3
"""
parse_popcorn.py - collect Popcorn (popcorn fit) results into the gencorr
standard results table.

Inputs (one cross-ancestry pair per Popcorn run, as written by POPCORN_FIT):
  --results  <pair>.popcorn          Popcorn's table: rows h1^2, h2^2 and pgi
                                     (genetic impact) or pge (genetic effect);
                                     columns 'Val (obs)', SE, ['Val (Lia)'], Z,
                                     'P (Z)'. Empty when Popcorn failed.
  --logs     <pair>.popcorn.o.log    Popcorn's log: 'Invoking command: ...
                                     --sfile1 <trait1>.popcorn.tsv --sfile2
                                     <trait2>.popcorn.tsv', SNP counts, version
             <pair>.popcorn.e.log    Popcorn's stderr (error message on failure)

Output (--out prefix): <prefix>.popcorn.tsv in gencorr_common.RESULT_COLUMNS.

IMPORTANT: Popcorn's own Z and p for the correlation test pg < 1 ("is the
cross-ancestry correlation lower than perfect?"), not pg != 0. They are kept in
extra as z_pg_lt_1 / p_pg_lt_1. The standard z / p columns test pg != 0, like
every other method: z = pg / SE, two-sided normal p (computed here).
rg, its SE and both heritabilities are Popcorn's own numbers.
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

SUMSTATS_SUFFIX = ".popcorn.tsv"


def trait_ids(log_text: str, log_name: str) -> tuple[str, str]:
    ids = []
    for flag in ("--sfile1", "--sfile2"):
        m = re.search(rf"{flag} (\S+)", log_text)
        if not m or not m.group(1).endswith(SUMSTATS_SUFFIX):
            die(f"{log_name}: no '{flag} <id>{SUMSTATS_SUFFIX}' in the logged command")
        ids.append(Path(m.group(1)).name[: -len(SUMSTATS_SUFFIX)])
    return ids[0], ids[1]


def read_table(path: Path) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    # round_trip: values exactly as Popcorn printed them (default parser can differ in the last digit)
    return pd.read_csv(path, sep="\t", index_col=0, float_precision="round_trip")


def last_error(text: str) -> str | None:
    """Last line of a Python traceback / error message in Popcorn's stderr."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    errors = [line for line in lines if re.match(r"^\w*(Error|Exception)\b", line)]
    return errors[-1] if errors else (lines[-1] if lines else None)


def parse_pair(result_path: Path, out_log: Path, err_log: Path | None) -> dict:
    text = out_log.read_text(errors="replace")
    t1, t2 = trait_ids(text, out_log.name)
    tab = read_table(result_path)
    m = re.search(r"(\d+) SNPs remaining after filtering", text)
    n_snps = float(m.group(1)) if m else np.nan
    m = re.search(r"Popcorn version (\S+)", text)
    version = m.group(1) if m else None
    command = re.search(r"Invoking command: (.*)", text)
    command = command.group(1) if command else ""

    def get(row, col):
        try:
            return float(tab.loc[row, col])
        except (KeyError, TypeError, ValueError, AttributeError):
            return np.nan

    corr_row = None
    if tab is not None:
        corr_row = next((r for r in ("pge", "pgi", "pg") if r in tab.index), None)
    rg, rg_se = get(corr_row, "Val (obs)"), get(corr_row, "SE")

    if tab is None:
        status = "popcorn_error"
    elif corr_row is None or np.isnan(rg) or np.isnan(rg_se) or rg_se <= 0:
        status = "rg_not_estimable"
    else:
        status = "ok"
    z = rg / rg_se if status == "ok" else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if status == "ok" else np.nan

    extra = {
        "correlation": {"pge": "genetic_effect", "pgi": "genetic_impact"}.get(corr_row, corr_row),
        "fit": "mle" if "--use_mle" in command else "regression",
        "z_pg_lt_1": get(corr_row, "Z"), "p_pg_lt_1": get(corr_row, "P (Z)"),
        "p_method": "z = rg / SE, two-sided normal (tests rg != 0); Popcorn's own Z/P test rg < 1 (z_pg_lt_1, p_pg_lt_1)",
        "h2_1_liab": get("h1^2", "Val (Lia)"), "h2_2_liab": get("h2^2", "Val (Lia)"),
        "popcorn_version": version,
        "error": last_error(err_log.read_text(errors="replace")) if (status == "popcorn_error" and err_log and err_log.exists()) else None,
        "log": out_log.name,
    }
    extra = {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in extra.items()}

    return {
        "method": "popcorn", "scope": "global", "trait1": t1, "trait2": t2,
        "locus": np.nan, "chr": np.nan, "start": np.nan, "stop": np.nan, "n_snps": n_snps,
        "rg": rg, "rg_se": rg_se,
        "rg_ci_low": rg - 1.96 * rg_se if status == "ok" else np.nan,
        "rg_ci_high": rg + 1.96 * rg_se if status == "ok" else np.nan,
        "z": z, "p": p,
        "h2_1": get("h1^2", "Val (obs)"), "h2_1_se": get("h1^2", "SE"),
        "h2_2": get("h2^2", "Val (obs)"), "h2_2_se": get("h2^2", "SE"),
        "status": status, "extra": json.dumps(extra, sort_keys=True),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", nargs="*", default=[], help="<pair>.popcorn files")
    p.add_argument("--logs", nargs="+", required=True, help="<pair>.popcorn.o.log and .e.log files")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    results = {Path(r).name: Path(r) for r in args.results}
    logs = {Path(lg).name: Path(lg) for lg in args.logs}
    rows = []
    for name in sorted(n for n in logs if n.endswith(".popcorn.o.log")):
        stem = name[: -len(".o.log")]                      # <pair>.popcorn
        rows.append(parse_pair(results.get(stem, Path("/nonexistent")), logs[name],
                               logs.get(f"{stem}.e.log")))
    if not rows:
        die("no <pair>.popcorn.o.log files given")
    df = validate_results(pd.DataFrame(rows, columns=RESULT_COLUMNS), "popcorn")
    write_tsv(df, f"{args.out}.popcorn.tsv")
    print(f"popcorn: {len(df)} pairs, {int((df['status'] == 'ok').sum())} ok")


if __name__ == "__main__":
    main()
