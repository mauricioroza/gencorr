#!/usr/bin/env python3
"""
ldsc_coverage.py - how much of the LD score reference one trait covers.

Genome-wide estimators (LDSC, SumHer) assume summary statistics for SNPs
spread over the whole genome. A targeted array, an exome-imputed GWAS, a
regional file or one with missing chromosomes can still give an h2 and rg,
but they are not meaningful. This script measures coverage against the
reference the regression actually uses, so the confidence tier
(gencorr_common.row_confidence) does not depend on the size of that
reference:

  snp_fraction    SNPs with a Z-score / reference SNPs. The reference is the
                  LD score SNP set, restricted to the --merge-alleles list when
                  munge used one (munge then writes every list SNP, with an
                  empty Z when the trait lacks it).
  block_coverage  share of reference blocks in which the trait has at least
                  BLOCK_MIN_FRACTION of the block's SNPs. Blocks are runs of
                  ~BLOCK_SIZE consecutive reference SNPs (by position) on each
                  chromosome, so they follow the reference's own SNP density
                  and need no build-specific region file.

Output: a small "key: value" text file read by parse_ldsc.py.

Usage:
  ldsc_coverage.py --sumstats <id>.sumstats.gz --ref-ld ref_ld/<prefix>
                   [--split-by-chr] --out <id>.coverage.log
"""

import argparse
import glob
import sys

import numpy as np
import pandas as pd

BLOCK_SIZE = 1000          # reference SNPs per block (~2.5 Mb with HapMap3)
BLOCK_MIN_FRACTION = 0.25  # a block counts as covered from this share of its SNPs


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def ld_files(prefix: str, split_by_chr: bool) -> list[str]:
    pattern = f"{prefix}[0-9]*.l2.ldscore*" if split_by_chr else f"{prefix}.l2.ldscore*"
    files = sorted(glob.glob(pattern))
    if not files:
        die(f"no LD score files match {pattern}")
    return files


def read_reference(prefix: str, split_by_chr: bool) -> pd.DataFrame:
    parts = [pd.read_csv(f, sep=r"\s+", usecols=["CHR", "SNP", "BP"],
                         dtype={"CHR": str, "SNP": str})
             for f in ld_files(prefix, split_by_chr)]
    ref = pd.concat(parts, ignore_index=True).drop_duplicates("SNP")
    return ref


def assign_blocks(ref: pd.DataFrame, block_size: int = BLOCK_SIZE) -> pd.Series:
    """Block label per reference SNP: each chromosome is cut, by position, into
    round(n / block_size) (at least 1) blocks of equal SNP count."""
    ref = ref.sort_values(["CHR", "BP"], kind="mergesort")
    labels = pd.Series(index=ref.index, dtype=object)
    for chrom, idx in ref.groupby("CHR", sort=False).groups.items():
        n = len(idx)
        n_blocks = max(1, int(round(n / block_size)))
        rank = np.arange(n)
        labels.loc[idx] = [f"{chrom}:{b}" for b in (rank * n_blocks) // n]
    return labels


def coverage(ref: pd.DataFrame, sumstats: pd.DataFrame,
             block_size: int = BLOCK_SIZE, min_fraction: float = BLOCK_MIN_FRACTION) -> dict:
    listed = sumstats["SNP"]
    used = set(sumstats.loc[sumstats["Z"].notna(), "SNP"])
    if sumstats["Z"].isna().any():          # munge --merge-alleles: the list is the universe
        ref = ref[ref["SNP"].isin(set(listed))]
    if ref.empty:
        die("no reference SNP left to measure coverage against")
    ref = ref.assign(block=assign_blocks(ref, block_size), used=ref["SNP"].isin(used))
    per_block = ref.groupby("block")["used"].mean()
    return {
        "reference_snps": int(len(ref)),
        "snps_used": int(ref["used"].sum()),
        "snp_fraction": float(ref["used"].mean()),
        "blocks": int(len(per_block)),
        "block_coverage": float((per_block >= min_fraction).mean()),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sumstats", required=True, help="munged <id>.sumstats.gz (SNP, Z columns)")
    p.add_argument("--ref-ld", required=True, help="LD score prefix, as given to ldsc.py")
    p.add_argument("--split-by-chr", action="store_true")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    sumstats = pd.read_csv(args.sumstats, sep="\t", usecols=["SNP", "Z"], dtype={"SNP": str})
    cov = coverage(read_reference(args.ref_ld, args.split_by_chr), sumstats)
    with open(args.out, "w") as fh:
        fh.write("gencorr ldsc_coverage.py: coverage of the LD score reference\n")
        fh.write(f"Block size: {BLOCK_SIZE} reference SNPs; covered from {BLOCK_MIN_FRACTION:g} of them\n")
        for key, value in cov.items():
            fh.write(f"{key}: {value:.6g}\n" if isinstance(value, float) else f"{key}: {value}\n")
    print(open(args.out).read(), end="")


if __name__ == "__main__":
    main()
