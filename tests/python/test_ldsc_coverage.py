"""
Tests for bin/ldsc_coverage.py (share of the LD score reference a trait covers)
and for how parse_ldsc.py / parse_sumher.py pass coverage to the confidence tier.
"""

import gzip
import json

import numpy as np
import pandas as pd

from conftest import BIN, run

import ldsc_coverage as lc
import parse_ldsc
import parse_sumher


def reference(n_per_chr=4000, chroms=("1", "2")):
    return pd.DataFrame({
        "CHR": [c for c in chroms for _ in range(n_per_chr)],
        "SNP": [f"rs{c}_{i}" for c in chroms for i in range(n_per_chr)],
        "BP": [1000 * (i + 1) for _ in chroms for i in range(n_per_chr)],
    })


def sumstats(snps, listed=None):
    """Munged file: SNPs in `snps` have Z; others in `listed` have an empty Z
    (what munge --merge-alleles writes)."""
    listed = list(listed) if listed is not None else list(snps)
    have = set(snps)
    return pd.DataFrame({"SNP": listed, "Z": [1.0 if s in have else np.nan for s in listed]})


def test_full_coverage():
    ref = reference()
    cov = lc.coverage(ref, sumstats(ref["SNP"]))
    assert cov["snp_fraction"] == 1 and cov["block_coverage"] == 1
    assert cov["blocks"] == 8 and cov["reference_snps"] == 8000      # 4 blocks of 1,000 per chr


def test_even_thinning_keeps_block_coverage():
    # every third SNP: fraction 1/3 everywhere, every block still >= 25%
    ref = reference()
    cov = lc.coverage(ref, sumstats(ref["SNP"].iloc[::3]))
    assert abs(cov["snp_fraction"] - 1 / 3) < 0.01 and cov["block_coverage"] == 1


def test_missing_chromosome_lowers_block_coverage():
    ref = reference()
    cov = lc.coverage(ref, sumstats(ref.loc[ref.CHR == "1", "SNP"]))
    assert cov["snp_fraction"] == 0.5 and cov["block_coverage"] == 0.5


def test_clustered_snps_lower_block_coverage():
    # exome-like: all SNPs of a few blocks, nothing elsewhere
    ref = reference()
    cov = lc.coverage(ref, sumstats(ref["SNP"].iloc[:1500]))
    assert cov["block_coverage"] == 2 / 8


def test_merge_alleles_list_is_the_denominator():
    # reference has 8,000 SNPs, the merge list only chr 1; the trait has all of chr 1
    ref = reference()
    chr1 = ref.loc[ref.CHR == "1", "SNP"]
    listed = list(chr1) + ["rs_not_in_ref"]
    cov = lc.coverage(ref, sumstats(chr1, listed=listed + ["rs_missing"]))
    assert cov["reference_snps"] == 4000 and cov["snp_fraction"] == 1


def test_cli_split_by_chr(tmp_path):
    ref = reference()
    (tmp_path / "ref").mkdir()
    for c, part in ref.groupby("CHR"):
        part.assign(L2=1.0)[["CHR", "SNP", "BP", "L2"]].to_csv(
            tmp_path / "ref" / f"{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
    with gzip.open(tmp_path / "t.sumstats.gz", "wt") as fh:
        sumstats(ref["SNP"].iloc[::2]).assign(A1="A", A2="G", N=100).to_csv(fh, sep="\t", index=False)
    run(["python3", BIN / "ldsc_coverage.py", "--sumstats", tmp_path / "t.sumstats.gz",
         "--ref-ld", f"{tmp_path}/ref/", "--split-by-chr", "--out", tmp_path / "t.coverage.log"])
    parsed = parse_ldsc.parse_coverage_log(tmp_path / "t.coverage.log")
    assert parsed == {"snp_fraction": 0.5, "block_coverage": 1.0}


def test_parse_ldsc_puts_coverage_in_h2_table_and_extra(tmp_path):
    (tmp_path / "a.coverage.log").write_text("snp_fraction: 0.13\nblock_coverage: 0.2\n")
    h2 = parse_ldsc.collect_h2([tmp_path / "a.coverage.log"])
    row = h2.set_index("trait").loc["a"]
    assert row["snp_fraction"] == 0.13 and row["block_coverage"] == 0.2


def test_sumher_fraction_of_tagging_predictors():
    log = ("Reading details for 90418 predictors from small.tagging\n"
           "Warning, bmi.sumher.tsv contains summary statistics for only 39310 of the 90418 predictors\n")
    assert parse_sumher.tagging_predictors(log) == 90418
    assert np.isnan(parse_sumher.tagging_predictors("nothing here"))
