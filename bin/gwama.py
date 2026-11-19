#!/usr/bin/env python3
"""
gwama.py - N-weighted multivariate GWAMA of genetically correlated traits.

Combines the harmonized summary statistics of K traits into ONE new trait, in
the same harmonized format, so it can be analysed by every method like any
other trait. Method: Baselmans et al. (2019) Nat Genet 51:445-451, N-weighted
GWAMA ("N-GWAMA"), reimplemented from the published formula and the authors'
R function (github.com/baselmans/multivariate_GWAMA, version 1_2_6); the oracle
test in tests/python/test_gwama.py compares the two.

Per SNP, with w_i = N_i * h2_i for the traits that have the SNP:

    Z  = sum_i sqrt(w_i) Z_i / sqrt( sum_ij sqrt(w_i) sqrt(w_j) C_ij )
    N  = sqrt(N)' C^-1 sqrt(N)                      ("effective" N, the default)
    or   sum_i N_i                                  ("observed" N)
    EAF  = N-weighted mean of the traits' EAF
    BETA = Z / sqrt(N_obs 2 EAF (1 - EAF)),  SE = 1 / sqrt(N_obs 2 EAF (1 - EAF))
    P    = two-sided, from Z

C is the cross-trait intercept (CTI) matrix: univariate LDSC intercepts on the
diagonal, cross-trait LDSC intercepts (gcov_int) off the diagonal. It carries
the error correlation from sample overlap and the inflation from population
stratification. h2 is the observed-scale SNP heritability from LDSC.

Inputs
  --traits    trait IDs, in the order of --sumstats (the first trait that
              has a SNP decides its effect allele and position)
  --sumstats  harmonized files (bin/gencorr_common.py HARMONIZED_COLUMNS)
  --signs     +1 / -1 per trait (default all +1). -1 reverses a trait's
              effects before combining, e.g. neuroticism with well-being.
  either      --ldsc-h2 <gencorr.ldsc_h2.tsv> --ldsc-results <gencorr.ldsc.tsv>
              (the pipeline: h2, intercepts and gcov_int from LDSC)
  or          --cti <matrix file> --h2 <values>   (the user's own; used by
              the oracle test with the authors' example inputs)

Outputs (--out prefix)
  <prefix>.harmonized.tsv.gz   the combined trait (INFO is NA)
  <prefix>.gwama.json          inputs used (signs, h2, CTI, rg) and QC counts

SNPs are combined over the union of the traits (the R function's default
method='union'): a SNP missing in some traits is combined from the others.
Alleles are aligned to the reference trait; a SNP whose alleles do not match
(not even swapped) is dropped for that trait only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gencorr_common import HARMONIZED_COLUMNS, NA_VALUES, die, read_tsv  # noqa: E402

GENOME_WIDE = 5e-8


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def read_harmonized(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=NA_VALUES, keep_default_na=False,
                     dtype={"SNP": str, "CHR": str, "A1": str, "A2": str},
                     usecols=["SNP", "CHR", "BP", "A1", "A2", "Z", "N", "EAF"])
    return df


def read_cti_file(path: str | Path, traits: list[str]) -> np.ndarray:
    """A square matrix with trait names as row and column names (R write.table /
    read.table format: the header may or may not have a corner cell)."""
    mat = pd.read_csv(path, sep=r"\s+", index_col=0, engine="python")
    if list(mat.columns) != list(mat.index):
        die(f"{path}: row and column names of the CTI matrix differ")
    if len(mat) != len(traits):
        die(f"{path}: CTI matrix has {len(mat)} rows for {len(traits)} traits")
    # Matched by position when the names are not the trait IDs (the authors'
    # example names them T1..T4), by name otherwise.
    if set(mat.index) == set(traits):
        mat = mat.loc[traits, traits]
    return mat.to_numpy(dtype=float)


def cti_from_ldsc(h2_path: str | Path, results_path: str | Path, traits: list[str]
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (h2, CTI, rg) for the traits from parse_ldsc.py's tables.

    CTI diagonal = univariate LDSC intercept; off-diagonal = cross-trait
    intercept (gcov_int, averaged over both orientations when both were run).
    """
    h2 = read_tsv(h2_path)
    h2["trait"] = h2["trait"].astype(str)
    h2 = h2.set_index("trait")
    missing = [t for t in traits if t not in h2.index]
    if missing:
        die(f"no LDSC h2 for {missing} in {h2_path}")
    k = len(traits)
    h2_obs = np.array([h2.at[t, "h2_obs"] for t in traits], dtype=float)
    cti = np.full((k, k), np.nan)
    for i, t in enumerate(traits):
        cti[i, i] = h2.at[t, "intercept"]
        if not np.isfinite(cti[i, i]) or cti[i, i] <= 0:
            die(f"{t}: no usable LDSC h2 intercept ({h2.at[t, 'intercept']}); "
                "N-GWAMA needs the estimated (not constrained) intercept")

    res = read_tsv(results_path)
    res = res[res["method"] == "ldsc"]
    rg = np.eye(k)
    found: dict[tuple[int, int], list[tuple[float, float]]] = {}
    index = {t: i for i, t in enumerate(traits)}
    for _, r in res.iterrows():
        a, b = str(r["trait1"]), str(r["trait2"])
        if a not in index or b not in index or a == b:
            continue
        extra = json.loads(r["extra"]) if isinstance(r["extra"], str) else {}
        gcov = extra.get("gcov_int")
        gcov = np.nan if gcov is None else float(gcov)
        key = tuple(sorted((index[a], index[b])))
        found.setdefault(key, []).append((gcov, float(r["rg"])))
    for i in range(k):
        for j in range(i + 1, k):
            values = found.get((i, j))
            gcovs = [g for g, _ in values or [] if np.isfinite(g)]
            if not gcovs:
                die(f"no LDSC cross-trait intercept (gcov_int) for {traits[i]} x {traits[j]} "
                    f"in {results_path}")
            cti[i, j] = cti[j, i] = float(np.mean(gcovs))
            rgs = [x for _, x in values if np.isfinite(x)]
            rg[i, j] = rg[j, i] = float(np.mean(rgs)) if rgs else np.nan
    return h2_obs, cti, rg


# --------------------------------------------------------------------------- #
# Alignment
# --------------------------------------------------------------------------- #
def align(tables: list[pd.DataFrame]) -> tuple[pd.DataFrame, np.ndarray, np.ndarray,
                                               np.ndarray, list[dict]]:
    """Union of SNPs, alleles aligned to the first trait that has each SNP.

    Returns (snp_info, Z, N, EAF, per-trait counts); the matrices are
    n_snps x n_traits with NaN where a trait lacks the SNP.
    """
    counts = []
    clean = []
    for df in tables:
        n_in = len(df)
        dup = df["SNP"].duplicated(keep=False)
        df = df[~dup & df["Z"].notna() & df["N"].notna()].copy()
        df["A1"] = df["A1"].str.upper()
        df["A2"] = df["A2"].str.upper()
        clean.append(df)
        counts.append({"rows_in": n_in, "duplicated_or_missing": int(n_in - len(df))})

    ref = (pd.concat([d[["SNP", "CHR", "BP", "A1", "A2"]] for d in clean], ignore_index=True)
           .drop_duplicates("SNP", keep="first")
           .reset_index(drop=True))
    ref_index = pd.Index(ref["SNP"])
    n_snps, k = len(ref), len(clean)
    Z = np.full((n_snps, k), np.nan)
    N = np.full((n_snps, k), np.nan)
    EAF = np.full((n_snps, k), np.nan)
    for j, df in enumerate(clean):
        rows = ref_index.get_indexer(df["SNP"])
        ra1 = ref["A1"].to_numpy()[rows]
        ra2 = ref["A2"].to_numpy()[rows]
        a1 = df["A1"].to_numpy()
        a2 = df["A2"].to_numpy()
        same = (a1 == ra1) & (a2 == ra2)
        swap = (a1 == ra2) & (a2 == ra1)
        keep = same | swap
        sign = np.where(swap, -1.0, 1.0)
        z = df["Z"].to_numpy(dtype=float) * sign
        eaf = df["EAF"].to_numpy(dtype=float)
        eaf = np.where(swap, 1.0 - eaf, eaf)
        Z[rows[keep], j] = z[keep]
        N[rows[keep], j] = df["N"].to_numpy(dtype=float)[keep]
        EAF[rows[keep], j] = eaf[keep]
        counts[j].update({"aligned": int(same.sum()), "flipped": int((swap & ~same).sum()),
                          "allele_mismatch": int((~keep).sum())})
    return ref, Z, N, EAF, counts


# --------------------------------------------------------------------------- #
# N-GWAMA
# --------------------------------------------------------------------------- #
def _patterns(present: np.ndarray):
    """Group rows by which traits are present: yields (row mask, trait index)."""
    k = present.shape[1]
    code = present.astype(np.int64) @ (1 << np.arange(k, dtype=np.int64))
    for c in np.unique(code):
        idx = np.array([(c >> j) & 1 for j in range(k)], dtype=bool)
        if idx.any():
            yield code == c, idx


def effective_n(N: np.ndarray, cti: np.ndarray) -> np.ndarray:
    """sqrt(N)' C^-1 sqrt(N) per SNP, over the traits that have it. As in the R
    function: plain sum when a SNP has one trait or no cross-trait intercepts
    (off-diagonal sum 0, then the diagonal intercepts are not used either)."""
    present = ~np.isnan(N)
    out = np.full(N.shape[0], np.nan)
    lower = np.tril_indices(cti.shape[0], -1)
    if cti[lower].sum() == 0:
        return np.nansum(N, axis=1)
    for rows, idx in _patterns(present):
        sub = cti[np.ix_(idx, idx)]
        n = N[np.ix_(rows, idx)]
        if idx.sum() == 1:
            out[rows] = n[:, 0]
        elif sub[np.tril_indices(idx.sum(), -1)].sum() == 0:
            out[rows] = n.sum(axis=1)
        else:
            s = np.sqrt(n)
            out[rows] = np.einsum("ri,ij,rj->r", s, np.linalg.inv(sub), s)
    return out


def ngwama(Z: np.ndarray, N: np.ndarray, EAF: np.ndarray, h2: np.ndarray, cti: np.ndarray,
           n_out: str = "effective") -> pd.DataFrame:
    sqrt_w = np.sqrt(N * h2[np.newaxis, :])             # NaN where the trait lacks the SNP
    s = np.nan_to_num(sqrt_w)
    numerator = np.nansum(Z * sqrt_w, axis=1)
    variance = np.einsum("ri,ij,rj->r", s, cti, s)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = numerator / np.sqrt(variance)
    p = 2 * stats.norm.sf(np.abs(z))
    n_obs = np.nansum(N, axis=1)
    has_eaf = ~np.isnan(EAF)
    n_eaf = np.where(has_eaf, N, 0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        eaf = np.where(has_eaf, EAF * N, 0).sum(axis=1) / n_eaf
        eaf = np.where(n_eaf > 0, eaf, np.nan)
        se = 1 / np.sqrt(n_obs * 2 * eaf * (1 - eaf))
        beta = z * se
    se = np.where(np.isfinite(se), se, np.nan)
    beta = np.where(np.isfinite(beta), beta, np.nan)
    n = effective_n(N, cti) if n_out == "effective" else n_obs
    return pd.DataFrame({"Z": z, "BETA": beta, "SE": se, "P": p, "N": n, "N_obs": n_obs,
                         "EAF": eaf})


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def check_inputs(traits: list[str], h2: np.ndarray, cti: np.ndarray, rg: np.ndarray | None,
                 min_rg: float) -> list[str]:
    warnings = []
    for t, v in zip(traits, h2):
        if not np.isfinite(v) or v <= 0:
            die(f"{t}: SNP heritability {v} is not positive; N-GWAMA weights traits by N * h2")
        if v > 1:
            warnings.append(f"{t}: SNP heritability {v} is larger than 1")
    if not np.allclose(cti, cti.T):
        die("the CTI matrix is not symmetric")
    min_eig = float(np.linalg.eigvalsh(cti).min())
    if min_eig <= 0:
        die(f"the CTI matrix is not positive definite (smallest eigenvalue {min_eig:.3g}); "
            "check the LDSC intercepts of the group's traits")
    if rg is not None:
        k = len(traits)
        for i in range(k):
            for j in range(i + 1, k):
                r = rg[i, j]
                pair = f"{traits[i]} x {traits[j]}"
                if not np.isfinite(r):
                    warnings.append(f"{pair}: LDSC did not estimate rg")
                elif r < 0:
                    warnings.append(f"{pair}: rg = {r:.3f} is negative after signs; "
                                    "give one of them sign -1 in the groups file")
                elif r < min_rg:
                    warnings.append(f"{pair}: rg = {r:.3f} < {min_rg}; Baselmans et al. found "
                                    "no gain over single-trait analysis at rg <= 0.4")
    return warnings


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_signs(values: list[str] | None, k: int) -> np.ndarray:
    if not values:
        return np.ones(k)
    if len(values) != k:
        die(f"--signs has {len(values)} values for {k} traits")
    table = {"1": 1.0, "+1": 1.0, "+": 1.0, "-1": -1.0, "-": -1.0}
    bad = [v for v in values if v not in table]
    if bad:
        die(f"--signs: values must be 1 or -1 (got {bad})")
    return np.array([table[v] for v in values])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--id", required=True, help="ID of the combined trait")
    p.add_argument("--traits", nargs="+", required=True)
    p.add_argument("--sumstats", nargs="+", required=True)
    p.add_argument("--signs", nargs="+")
    p.add_argument("--ldsc-h2")
    p.add_argument("--ldsc-results")
    p.add_argument("--cti")
    p.add_argument("--h2", nargs="+", type=float)
    p.add_argument("--n-out", choices=["effective", "observed"], default="effective",
                   help="N column of the output (default: effective)")
    p.add_argument("--min-rg", type=float, default=0.5,
                   help="warn when a pair's LDSC rg (after signs) is below this")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    traits, k = args.traits, len(args.traits)
    if k < 2:
        die("N-GWAMA needs at least 2 traits")
    if len(set(traits)) != k:
        die(f"--traits has duplicates: {traits}")
    if len(args.sumstats) != k:
        die(f"--sumstats has {len(args.sumstats)} files for {k} traits")
    signs = parse_signs(args.signs, k)

    if args.cti:
        if not args.h2 or len(args.h2) != k:
            die("--cti needs --h2 with one value per trait")
        h2, cti, rg = np.array(args.h2, dtype=float), read_cti_file(args.cti, traits), None
    elif args.ldsc_h2 and args.ldsc_results:
        h2, cti, rg = cti_from_ldsc(args.ldsc_h2, args.ldsc_results, traits)
    else:
        die("give either --ldsc-h2 and --ldsc-results, or --cti and --h2")

    # Reversing a trait's effects reverses its covariances with the others.
    cti = cti * np.outer(signs, signs)
    if rg is not None:
        rg = rg * np.outer(signs, signs)
    warnings = check_inputs(traits, h2, cti, rg, args.min_rg)
    for w in warnings:
        print(f"WARNING ({args.id}): {w}", file=sys.stderr)

    tables = []
    for t, path, s in zip(traits, args.sumstats, signs):
        df = read_harmonized(path)
        df["Z"] = df["Z"] * s
        tables.append(df)
    snp_info, Z, N, EAF, counts = align(tables)
    combined = ngwama(Z, N, EAF, h2, cti, args.n_out)

    out = pd.DataFrame({
        "SNP": snp_info["SNP"], "CHR": snp_info["CHR"], "BP": snp_info["BP"].astype("Int64"),
        "A1": snp_info["A1"], "A2": snp_info["A2"],
        "Z": combined["Z"], "BETA": combined["BETA"], "SE": combined["SE"], "P": combined["P"],
        "N": combined["N"], "EAF": combined["EAF"], "INFO": np.nan,
    })[HARMONIZED_COLUMNS]
    ok = out["Z"].notna() & np.isfinite(out["N"])
    out = out[ok]
    if out.empty:
        die(f"{args.id}: no SNPs left after combining")
    out.to_csv(f"{args.out}.harmonized.tsv.gz", sep="\t", index=False, na_rep="NA",
               float_format="%.8g", compression={"method": "gzip", "mtime": 0})

    present = ~np.isnan(Z[ok.to_numpy()])
    qc = {
        "id": args.id,
        "method": "N-weighted GWAMA (Baselmans et al. 2019)",
        "traits": [{"trait": t, "sign": int(s), "h2": float(v), **c}
                   for t, s, v, c in zip(traits, signs, h2, counts)],
        "cti": {"traits": traits, "matrix": np.round(cti, 6).tolist(),
                "source": "user" if args.cti else "ldsc"},
        "rg": None if rg is None else {"traits": traits, "matrix": np.round(rg, 6).tolist()},
        "n_out": args.n_out,
        "rows_out": int(len(out)),
        "snps_in_all_traits": int(present.all(axis=1).sum()),
        "snps_by_n_traits": {str(i): int((present.sum(axis=1) == i).sum()) for i in range(1, k + 1)},
        "genome_wide_significant": int((out["P"] <= GENOME_WIDE).sum()),
        "warnings": warnings,
    }
    with open(f"{args.out}.gwama.json", "w") as fh:
        json.dump(qc, fh, indent=2)
    print(json.dumps({k_: v for k_, v in qc.items() if k_ not in ("cti", "rg")}, indent=2))


if __name__ == "__main__":
    main()
