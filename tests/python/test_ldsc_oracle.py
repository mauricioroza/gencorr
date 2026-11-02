"""
Oracle tests: gencorr's LDSC path vs running LDSC 1.0.1 directly.

These answer "does the pipeline give the same numbers as the tool?":

1. Raw file -> munge_sumstats.py -> ldsc --rg     (what BigaGWAS / a user does)
   Raw file -> harmonize -> export -> munge -> rg (what gencorr does)
   must give identical rg, se, z, p and intercepts.

2. Batching pairs into one `ldsc --rg a,b,c` call (gencorr's LDSC_RG) must give
   the same rows as one call per pair.
"""

import pandas as pd
import pytest

from conftest import BIN, requires_ldsc, run

import parse_ldsc

pytestmark = requires_ldsc


def rg_table(log):
    t = parse_ldsc.parse_rg_log(log)
    return t[["trait1", "trait2", "rg", "se", "z", "p", "h2_obs", "h2_int", "gcov_int", "n_snps"]]


def ld_args(test_ld):
    return ["--ref-ld", test_ld / "all", "--w-ld", test_ld / "all"]


@pytest.fixture(scope="module")
def munged(vignette, test_ld, tmp_path_factory):
    """Both routes for bmi and neuro, munged with the same merge-alleles list."""
    d = tmp_path_factory.mktemp("ldsc")
    merge = ["--merge-alleles", test_ld / "merge_alleles.snplist.gz"]
    (d / "direct").mkdir()
    (d / "gencorr").mkdir()
    # Direct: the column choices a user would make by hand for these files.
    # MAF/INFO_UKB/STAT are ignored because the harmonized route does not carry
    # them either (BigaGWAS likewise passed --ignore for eaf).
    run(["munge_sumstats.py", "--sumstats", vignette / "bmi.sumstats.txt", "--out", d / "direct/bmi",
         "--snp", "SNPID_UKB", "--N-col", "NMISS", "--signed-sumstats", "BETA,0",
         "--ignore", "MAF,STAT,SE,INFO_UKB,CHR,BP", *merge])
    run(["munge_sumstats.py", "--sumstats", vignette / "neuro.sumstats.txt", "--out", d / "direct/neuro",
         "--N-col", "N_analyzed", "--signed-sumstats", "Z,0", "--ignore", "BETA,SE,CHR,POS", *merge])
    for t in ("bmi", "neuro"):
        run(["python3", BIN / "harmonize_sumstats.py", "--input", vignette / f"{t}.sumstats.txt",
             "--id", t, "--out", d / f"gencorr/{t}"])
        run(["python3", BIN / "export_sumstats.py", "--input", d / f"gencorr/{t}.harmonized.tsv.gz",
             "--format", "ldsc", "--out", d / f"gencorr/{t}.ldsc.tsv"])
        run(["munge_sumstats.py", "--sumstats", d / f"gencorr/{t}.ldsc.tsv.gz", "--out", d / f"gencorr/{t}",
             "--signed-sumstats", "Z,0", "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
             "--N-col", "N", *merge])
    return d


def test_gencorr_route_equals_direct_ldsc(munged, test_ld):
    for route in ("direct", "gencorr"):
        run(["ldsc.py", "--rg", "bmi.sumstats.gz,neuro.sumstats.gz", *ld_args(test_ld),
             "--out", "bmi__x.rg"], cwd=munged / route)
    direct = rg_table(munged / "direct/bmi__x.rg.log")
    ours = rg_table(munged / "gencorr/bmi__x.rg.log")
    pd.testing.assert_frame_equal(direct.reset_index(drop=True), ours.reset_index(drop=True))
    assert direct.rg.notna().all()


def test_rg_batching_does_not_change_estimates(munged, test_ld, vignette):
    d = munged / "gencorr"
    # add a third trait so a batch has two partners
    run(["python3", BIN / "harmonize_sumstats.py", "--input", vignette / "depression.sumstats.txt",
         "--id", "depression", "--out", d / "depression"])
    run(["python3", BIN / "export_sumstats.py", "--input", d / "depression.harmonized.tsv.gz",
         "--format", "ldsc", "--out", d / "depression.ldsc.tsv"])
    run(["munge_sumstats.py", "--sumstats", d / "depression.ldsc.tsv.gz", "--out", d / "depression",
         "--signed-sumstats", "Z,0", "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
         "--N-col", "N", "--merge-alleles", test_ld / "merge_alleles.snplist.gz"])
    run(["ldsc.py", "--rg", "bmi.sumstats.gz,neuro.sumstats.gz,depression.sumstats.gz",
         *ld_args(test_ld), "--out", "bmi__batch.rg"], cwd=d)
    run(["ldsc.py", "--rg", "bmi.sumstats.gz,neuro.sumstats.gz", *ld_args(test_ld),
         "--out", "bmi__s1.rg"], cwd=d)
    run(["ldsc.py", "--rg", "bmi.sumstats.gz,depression.sumstats.gz", *ld_args(test_ld),
         "--out", "bmi__s2.rg"], cwd=d)
    batch = rg_table(d / "bmi__batch.rg.log").reset_index(drop=True)
    single = pd.concat([rg_table(d / "bmi__s1.rg.log"), rg_table(d / "bmi__s2.rg.log")]).reset_index(drop=True)
    pd.testing.assert_frame_equal(batch, single)
