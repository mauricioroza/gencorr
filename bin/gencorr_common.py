"""
Shared definitions for the gencorr helper scripts.

Everything that more than one script must agree on lives here:

* HARMONIZED_COLUMNS  - the column contract of a harmonized summary statistics
                        file (output of harmonize_sumstats.py, input of every
                        method's export step).
* RESULT_COLUMNS      - the column contract of the standard results table that
                        every method must produce. merge_results.py refuses
                        tables that do not follow it.

If you add a method, you only need RESULT_COLUMNS (and optionally an entry in
export_sumstats.FORMATS). See docs/extending.md.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import sys
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Harmonized summary statistics
# --------------------------------------------------------------------------- #
# A1 is ALWAYS the effect allele; Z, BETA and EAF refer to A1.
# CHR/BP/BETA/SE/EAF/INFO may be NA when the source file does not have them.
HARMONIZED_COLUMNS = [
    "SNP",   # variant ID as given in the source (usually rsID), lower-cased "rs"
    "CHR",   # chromosome, string without "chr" prefix ("1".."22","X")
    "BP",    # base-pair position (build declared in the samplesheet)
    "A1",    # effect allele, upper case, one of A/C/G/T
    "A2",    # other allele,  upper case, one of A/C/G/T
    "Z",     # signed Z-score for A1
    "BETA",  # effect size for A1 (log-odds for binary traits), may be NA
    "SE",    # standard error of BETA, may be NA
    "P",     # two-sided p-value, recomputed from Z when missing or zero
    "N",     # per-variant total sample size
    "EAF",   # frequency of A1, may be NA
    "INFO",  # imputation quality, may be NA
]

# --------------------------------------------------------------------------- #
# Standard results table
# --------------------------------------------------------------------------- #
# One row per (method, trait pair, locus). For global methods locus/chr/start/
# stop are NA. Columns a method cannot fill are NA; method-specific values go
# into `extra` as a JSON object so the table stays rectangular.
RESULT_COLUMNS = [
    "method",      # e.g. "ldsc", "lava"
    "scope",       # "global" or "local"
    "trait1",
    "trait2",
    "locus",       # locus ID (local methods only)
    "chr",
    "start",
    "stop",
    "n_snps",      # SNPs used for this estimate, when the tool reports it
    "rg",          # genetic correlation estimate
    "rg_se",       # standard error (NA for LAVA, which reports CIs)
    "rg_ci_low",
    "rg_ci_high",
    "z",
    "p",
    "h2_1",        # heritability of trait1 (global h2 or local h2, observed scale)
    "h2_1_se",
    "h2_2",
    "h2_2_se",
    "status",      # "ok" or a short reason the estimate is missing
    "extra",       # JSON object with method-specific fields
]

VALID_SCOPES = {"global", "local"}

NA_VALUES = ["NA", "na", "NaN", "nan", "", ".", "None", "null", "NULL"]


def open_text(path: str | Path):
    """Open plain or gzip-compressed text transparently (by magic bytes)."""
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def empty_results() -> pd.DataFrame:
    return pd.DataFrame(columns=RESULT_COLUMNS)


def validate_results(df: pd.DataFrame, source: str = "<table>") -> pd.DataFrame:
    """Check a results table against the contract; return it with ordered columns."""
    missing = [c for c in RESULT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{source}: missing result columns {missing}")
    extra_cols = [c for c in df.columns if c not in RESULT_COLUMNS]
    if extra_cols:
        raise ValueError(
            f"{source}: unexpected columns {extra_cols}; put method-specific "
            "values in the JSON 'extra' column"
        )
    bad_scope = set(df["scope"].dropna()) - VALID_SCOPES
    if bad_scope:
        raise ValueError(f"{source}: invalid scope values {sorted(bad_scope)}")
    if df["method"].isna().any() or df["trait1"].isna().any() or df["trait2"].isna().any():
        raise ValueError(f"{source}: method/trait1/trait2 must never be NA")
    for value in df["extra"].dropna():
        try:
            json.loads(value)
        except (TypeError, json.JSONDecodeError) as err:
            raise ValueError(f"{source}: 'extra' is not valid JSON: {value!r}") from err
    return df[RESULT_COLUMNS]


def write_tsv(df: pd.DataFrame, path: str | Path) -> None:
    """Write a TSV with 'NA' for missing values and no quoting.

    Quoting is disabled so the JSON 'extra' column is written as-is (JSON never
    contains tabs or newlines here). Read these files back with read_tsv().
    """
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.10g",
              quoting=csv.QUOTE_NONE, quotechar="\x1f")


def read_tsv(path: str | Path) -> pd.DataFrame:
    """Read a gencorr TSV (results or tables) written by write_tsv() or R."""
    return pd.read_csv(path, sep="\t", na_values=NA_VALUES, keep_default_na=False,
                       quoting=csv.QUOTE_NONE, dtype={"trait1": str, "trait2": str,
                                                      "locus": str, "chr": str})


def die(message: str, code: int = 1) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(code)
