#!/usr/bin/env python3
"""
harmonize_sumstats.py - turn one GWAS summary statistics file, in any of the
common layouts, into the gencorr harmonized format (see gencorr_common.py).

What it does, in order (every step is counted in the QC report):

 1. Detect the delimiter (tab, comma or whitespace) and read everything as text.
 2. Map source columns to standard names using a synonym table
    (COLUMN_SYNONYMS below); --colmap overrides any automatic choice.
 3. Normalise IDs and alleles: upper-case alleles, "rs" IDs lower-cased,
    "chr" prefix stripped from chromosomes.
 4. Derive BETA from OR (BETA = ln OR) when only OR is given.
 5. Choose the signed Z-score per variant, by priority (default --z-priority):
       z        an explicit Z column
       p        |Z| from the two-sided P, sign from BETA/OR
       beta_se  BETA / SE
    Using P first mirrors what LDSC munge_sumstats.py and LAVA do internally,
    so the Z passed to every method is the same Z the original tools would use.
 6. Resolve the per-variant N: N column > samplesheet --n > cases + controls.
 7. Drop variants that cannot be used: missing SNP/alleles/Z/N, non-SNV or
    multi-character alleles, A1 == A2, (optionally) strand-ambiguous A/T and C/G
    pairs, and every copy of a duplicated SNP ID.
 8. Optional QC thresholds (--min-info, --min-maf; off by default). Same rule
    as LDSC munge: a variant is removed when its value is <= the threshold.
    Variants with a missing INFO/MAF value are kept (genotyped SNPs often have
    no INFO). MAF comes from a MAF column if present, else min(EAF, 1 - EAF).
    If the file has neither column, the threshold cannot be applied and the
    QC report says so (no variant is removed for it).
 9. Recompute P from Z where P is missing or 0, and write the harmonized table
    plus a JSON QC report with counts and the median-chi2 lambda (lambda_GC).

It does NOT lift over genome builds or impute missing columns from a reference
panel. Without --min-info/--min-maf, thresholds are left to each method (LDSC
applies its own INFO > 0.9 / MAF > 0.01) so results stay comparable with
running the tools directly.
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
from gencorr_common import HARMONIZED_COLUMNS, NA_VALUES, die, open_text  # noqa: E402

# Synonyms are matched case-insensitively against the source header, in the
# order listed: the first header present wins. Keep the most specific names
# first. "STAT" is deliberately absent: depending on the GWAS software it is a
# Z, a t or a chi-square statistic, so map it explicitly (--colmap Z=STAT) if
# you know it is a Z.
COLUMN_SYNONYMS: dict[str, list[str]] = {
    "SNP": ["snp", "rsid", "rs_id", "hm_rsid", "variant_id", "markername",
            "snpid", "snpid_ukb", "rsid_ukb", "id", "marker"],
    "CHR": ["chr", "chrom", "chromosome", "hm_chrom", "#chrom", "#chr"],
    "BP": ["bp", "pos", "position", "base_pair_location", "hm_pos", "bp_hg19"],
    "A1": ["a1", "effect_allele", "hm_effect_allele", "ea", "alt", "allele1",
           "tested_allele"],
    "A2": ["a2", "other_allele", "hm_other_allele", "nea", "ref", "allele2",
           "non_effect_allele"],
    "Z": ["z", "zscore", "z_score", "z-score", "z_stat"],
    "BETA": ["beta", "hm_beta", "b", "effect", "logor", "log_or", "log_odds"],
    "OR": ["or", "odds_ratio", "hm_odds_ratio"],
    "SE": ["se", "standard_error", "stderr", "stderrlogor", "se_beta", "log_or_se"],
    "P": ["p", "pval", "p_value", "pvalue", "p-value", "p_val", "p.value"],
    "N": ["n", "n_total", "nmiss", "n_analyzed", "obs_ct", "n_samples",
          "sample_size", "neff", "n_eff"],
    "EAF": ["eaf", "effect_allele_frequency", "hm_effect_allele_frequency",
            "a1freq", "a1_freq", "freq", "frq", "freq1", "af"],
    "INFO": ["info", "info_ukb", "imputation_info", "impinfo", "r2"],
    # MAF is read only for --min-maf; it is not written to the harmonized file
    # (it has no allele orientation, unlike EAF)
    "MAF": ["maf", "minor_allele_frequency"],
}

EFFECTIVE_N_NAMES = {"neff", "n_eff"}
VALID_ALLELES = set("ACGT")
AMBIGUOUS_PAIRS = {frozenset("AT"), frozenset("CG")}
RS_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)
CHISQ_MEDIAN_1DF = stats.chi2.ppf(0.5, 1)  # 0.4549...
MIN_P = 1e-300  # floor so very strong signals are not lost to underflow


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="summary statistics (plain or .gz)")
    p.add_argument("--id", required=True, help="trait ID used for output names")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--n", type=float, help="total N if the file has no N column")
    p.add_argument("--n-cases", type=float, help="number of cases (binary traits)")
    p.add_argument("--n-controls", type=float, help="number of controls (binary traits)")
    p.add_argument("--colmap", default="",
                   help="comma-separated STANDARD=source overrides, e.g. 'Z=STAT,N=NMISS'")
    p.add_argument("--z-priority", default="z,p,beta_se",
                   help="order in which Z sources are tried (default: z,p,beta_se)")
    p.add_argument("--drop-ambiguous", action="store_true",
                   help="drop strand-ambiguous A/T and C/G SNPs")
    p.add_argument("--min-info", type=float, default=None,
                   help="remove variants with INFO <= this value (default: off)")
    p.add_argument("--min-maf", type=float, default=None,
                   help="remove variants with MAF <= this value (default: off)")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
# Reading and column mapping
# --------------------------------------------------------------------------- #
def detect_separator(path: Path) -> str:
    with open_text(path) as fh:
        header = fh.readline()
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    return r"\s+"


def read_table(path: Path) -> pd.DataFrame:
    sep = detect_separator(path)
    return pd.read_csv(path, sep=sep, dtype=str, na_values=NA_VALUES,
                       keep_default_na=False, engine="python" if sep == r"\s+" else "c",
                       compression="infer")


def map_columns(header: list[str], overrides: dict[str, str]) -> dict[str, str]:
    """Return {standard_name: source_column}."""
    lower_to_source = {}
    for col in header:
        lower_to_source.setdefault(col.strip().lower(), col)
    mapping: dict[str, str] = {}
    for standard, synonyms in COLUMN_SYNONYMS.items():
        if standard in overrides:
            continue
        for syn in synonyms:
            if syn in lower_to_source:
                mapping[standard] = lower_to_source[syn]
                break
    for standard, source in overrides.items():
        if standard not in COLUMN_SYNONYMS:
            die(f"--colmap: unknown standard column '{standard}' "
                f"(valid: {', '.join(COLUMN_SYNONYMS)})")
        if source not in header:
            die(f"--colmap: column '{source}' not found in header {header}")
        mapping[standard] = source
    # One source column must not feed two standard columns (e.g. 'b' as BETA and SE)
    used = {}
    for standard, source in mapping.items():
        if source in used:
            die(f"column '{source}' matched both {used[source]} and {standard}; "
                "use --colmap to disambiguate")
        used[source] = standard
    return mapping


def parse_overrides(text: str) -> dict[str, str]:
    out = {}
    for item in filter(None, (x.strip() for x in text.split(","))):
        if "=" not in item:
            die(f"--colmap entry '{item}' must look like STANDARD=source")
        key, value = item.split("=", 1)
        out[key.strip().upper()] = value.strip()
    return out


def to_num(series: pd.Series | None) -> pd.Series | None:
    """Numeric parse; always float so later NaN/float assignments are safe even
    when a column happens to contain only whole numbers (e.g. P = 0)."""
    if series is None:
        return None
    return pd.to_numeric(series, errors="coerce").astype(float)


# --------------------------------------------------------------------------- #
# Core harmonization (pure function, unit-tested in tests/python)
# --------------------------------------------------------------------------- #
def harmonize(raw: pd.DataFrame, mapping: dict[str, str], *, n=None, n_cases=None,
              n_controls=None, z_priority=("z", "p", "beta_se"),
              drop_ambiguous=False, min_info=None,
              min_maf=None) -> tuple[pd.DataFrame, dict]:
    qc: dict = {"rows_in": int(len(raw)), "column_mapping": mapping, "removed": {}}

    def col(name):
        return raw[mapping[name]] if name in mapping else None

    for required in ("SNP", "A1", "A2"):
        if required not in mapping:
            die(f"no {required} column found (header: {list(raw.columns)}); use --colmap")

    out = pd.DataFrame(index=raw.index)
    out["SNP"] = col("SNP").str.strip()
    is_rs = out["SNP"].str.match(RS_PATTERN, na=False)
    out.loc[is_rs, "SNP"] = out.loc[is_rs, "SNP"].str.lower()

    chrom = col("CHR")
    out["CHR"] = (chrom.str.strip().str.replace(r"^chr", "", regex=True, case=False)
                  if chrom is not None else np.nan)
    out["BP"] = to_num(col("BP")).round().astype("Int64") if "BP" in mapping else pd.NA

    out["A1"] = col("A1").str.strip().str.upper()
    out["A2"] = col("A2").str.strip().str.upper()

    beta = to_num(col("BETA"))
    odds = to_num(col("OR"))
    if beta is None and odds is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.log(odds.where(odds > 0))
        qc["beta_from_or"] = True
    se = to_num(col("SE"))
    pval = to_num(col("P"))
    zcol = to_num(col("Z"))

    # ---- choose Z per variant according to priority ---------------------- #
    z = pd.Series(np.nan, index=raw.index)
    source = pd.Series(pd.NA, index=raw.index, dtype="object")
    sign = np.sign(beta) if beta is not None else None
    for kind in z_priority:
        todo = z.isna()
        if kind == "z" and zcol is not None:
            cand = zcol
        elif kind == "p" and pval is not None and sign is not None:
            valid_p = pval.where((pval > 0) & (pval <= 1))
            cand = sign * stats.norm.isf(valid_p / 2)
            # a beta of exactly 0 gives sign 0 -> Z 0, which is what the data says
        elif kind == "beta_se" and beta is not None and se is not None:
            cand = beta / se.where(se > 0)
        elif kind not in {"z", "p", "beta_se"}:
            die(f"--z-priority: unknown source '{kind}'")
        else:
            continue
        cand = pd.Series(cand, index=raw.index).replace([np.inf, -np.inf], np.nan)
        take = todo & cand.notna()
        z[take] = cand[take]
        source[take] = kind
    qc["z_source_counts"] = {str(k): int(v) for k, v in source.value_counts().items()}
    out["Z"] = z
    out["BETA"] = beta if beta is not None else np.nan
    out["SE"] = se if se is not None else np.nan

    # ---- N --------------------------------------------------------------- #
    ncol = to_num(col("N"))
    if ncol is not None and mapping["N"].lower() in EFFECTIVE_N_NAMES:
        qc["warning_effective_n"] = (
            f"N taken from '{mapping['N']}' (effective sample size). LDSC "
            "liability-scale h2 assumes total N; check this is intended.")
    fallback_n = n
    if fallback_n is None and n_cases is not None and n_controls is not None:
        fallback_n = n_cases + n_controls
    if ncol is not None:
        out["N"] = ncol
        if fallback_n is not None:
            out["N"] = out["N"].fillna(fallback_n)
        qc["n_source"] = f"column:{mapping['N']}"
    elif fallback_n is not None:
        out["N"] = float(fallback_n)
        qc["n_source"] = "samplesheet"
    else:
        die("no N column and no n / n_cases+n_controls in the samplesheet")

    out["EAF"] = to_num(col("EAF")) if "EAF" in mapping else np.nan
    out["INFO"] = to_num(col("INFO")) if "INFO" in mapping else np.nan
    out["P"] = pval if pval is not None else np.nan

    # MAF only for --min-maf: MAF column, else derived from EAF; folded to <= 0.5
    maf = to_num(col("MAF")) if "MAF" in mapping else (
        out["EAF"] if "EAF" in mapping else None)
    out["_MAF"] = np.minimum(maf, 1 - maf) if maf is not None else np.nan
    qc["maf_source"] = (f"column:{mapping['MAF']}" if "MAF" in mapping else
                        f"EAF:{mapping['EAF']}" if "EAF" in mapping else None)

    # ---- filters, each counted ------------------------------------------- #
    def drop(mask: pd.Series, reason: str):
        nonlocal out
        mask = mask.fillna(False).astype(bool)
        qc["removed"][reason] = int(mask.sum())
        out = out.loc[~mask]

    drop(out["SNP"].isna() | (out["SNP"] == ""), "missing_snp")
    drop(~out["A1"].isin(VALID_ALLELES) | ~out["A2"].isin(VALID_ALLELES),
         "non_snv_or_missing_allele")
    drop(out["A1"] == out["A2"], "a1_equals_a2")
    drop(out["Z"].isna(), "missing_z")
    drop(out["N"].isna() | (out["N"] <= 0), "missing_or_nonpositive_n")
    if drop_ambiguous:
        pair = [frozenset(x) for x in zip(out["A1"], out["A2"])]
        drop(pd.Series([p in AMBIGUOUS_PAIRS for p in pair], index=out.index),
             "strand_ambiguous")
    # ---- optional QC thresholds (LDSC munge convention: remove value <= t) - #
    qc["thresholds"] = {"min_info": min_info, "min_maf": min_maf}
    if min_info is not None:
        if "INFO" in mapping:
            qc["info_missing_kept"] = int(out["INFO"].isna().sum())
            drop(out["INFO"] <= min_info, f"info_le_{min_info:g}")
        else:
            qc["warning_min_info"] = "--min-info requested but the file has no INFO column; not applied"
    if min_maf is not None:
        if qc["maf_source"] is not None:
            qc["maf_missing_kept"] = int(out["_MAF"].isna().sum())
            drop(out["_MAF"] <= min_maf, f"maf_le_{min_maf:g}")
        else:
            qc["warning_min_maf"] = "--min-maf requested but the file has no MAF or EAF column; not applied"

    dup = out["SNP"].duplicated(keep=False)
    drop(dup, "duplicated_snp_all_copies")

    # ---- P from Z where unusable ----------------------------------------- #
    need_p = out["P"].isna() | (out["P"] <= 0) | (out["P"] > 1)
    qc["p_recomputed_from_z"] = int(need_p.sum())
    out.loc[need_p, "P"] = 2 * stats.norm.sf(np.abs(out.loc[need_p, "Z"]))
    floored = out["P"] < MIN_P
    qc["p_floored_at_1e-300"] = int(floored.sum())
    out.loc[floored, "P"] = MIN_P

    chi2 = out["Z"] ** 2
    qc["rows_out"] = int(len(out))
    qc["lambda_gc"] = float(np.median(chi2) / CHISQ_MEDIAN_1DF) if len(out) else None
    qc["mean_chi2"] = float(chi2.mean()) if len(out) else None
    qc["median_n"] = float(out["N"].median()) if len(out) else None
    return out[HARMONIZED_COLUMNS].reset_index(drop=True), qc


def main(argv=None):
    args = parse_args(argv)
    raw = read_table(Path(args.input))
    mapping = map_columns(list(raw.columns), parse_overrides(args.colmap))
    harmonized, qc = harmonize(
        raw, mapping, n=args.n, n_cases=args.n_cases, n_controls=args.n_controls,
        z_priority=tuple(x.strip() for x in args.z_priority.split(",") if x.strip()),
        drop_ambiguous=args.drop_ambiguous, min_info=args.min_info,
        min_maf=args.min_maf)
    for key in ("warning_min_info", "warning_min_maf", "warning_effective_n"):
        if key in qc:
            print(f"WARNING ({args.id}): {qc[key]}", file=sys.stderr)
    qc = {"id": args.id, "input": Path(args.input).name, **qc}
    if qc["rows_out"] == 0:
        die(f"{args.id}: no variants left after harmonization; see QC: {json.dumps(qc)}")
    harmonized.to_csv(f"{args.out}.harmonized.tsv.gz", sep="\t", index=False,
                      na_rep="NA", float_format="%.8g",
                      compression={"method": "gzip", "mtime": 0})  # mtime 0: byte-identical reruns
    with open(f"{args.out}.qc.json", "w") as fh:
        json.dump(qc, fh, indent=2)
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
