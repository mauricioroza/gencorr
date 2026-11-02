"""
Oracle tests: gencorr's LAVA runner (bin/run_lava.R) vs LAVA's own API.

1. Same inputs, LAVA called as in its vignette (tests/oracle/lava_native.R)
   -> rho, CIs, p and local h2 must be IDENTICAL (same seeds, same draws).
2. Splitting loci into chunks must not change any number.
3. Harmonized inputs vs LAVA reading the raw vignette files itself
   -> same SNPs and statuses; neuro (Z column) local h2 identical;
      depression local h2 within 0.1% and rg within 1e-4.
   The one known difference: for variants with OR exactly 1 (LogOR 0; 435 in
   the depression file) LAVA's reader assigns a negative sign
   (ifelse(OR > 1, 1, -1)); the harmonizer gives them Z = 0 because there is
   no direction. LAVA projects the signed Z-scores onto LD principal
   components, so these signs move depression's local h2 (~1e-4 relative)
   and rg (5th decimal). The continuous trait is unaffected.
"""

import numpy as np
import pandas as pd
import pytest

from conftest import BIN, REPO, requires_lava, run

pytestmark = requires_lava

PAIR = ("neuro", "depression")        # continuous (Z column) x binary (LogOR/OR + P)
CASES = {"depression": (170756, 329443)}
SEED = 20201231
ORACLE = REPO / "tests" / "oracle" / "lava_native.R"


def write_info(path, files):
    rows = ["phenotype\tcases\tcontrols\tprevalence\tfilename"]
    for t, f in files.items():
        c, k = CASES.get(t, ("NA", "NA"))
        rows.append(f"{t}\t{c}\t{k}\tNA\t{f}")
    path.write_text("\n".join(rows) + "\n")


@pytest.fixture(scope="module")
def lava_inputs(vignette, tmp_path_factory):
    d = tmp_path_factory.mktemp("lava")
    files = {}
    for t in PAIR:
        run(["python3", BIN / "harmonize_sumstats.py", "--input", vignette / f"{t}.sumstats.txt",
             "--id", t, "--out", d / t])
        run(["python3", BIN / "export_sumstats.py", "--input", d / f"{t}.harmonized.tsv.gz",
             "--format", "lava", "--out", d / f"{t}.lava.tsv"])
        files[t] = d / f"{t}.lava.tsv"
    write_info(d / "info.txt", files)
    (d / "pairs.tsv").write_text(f"trait1\ttrait2\n{PAIR[0]}\t{PAIR[1]}\n")
    loci = (vignette / "test.loci").read_text().splitlines()
    (d / "loci.txt").write_text("\n".join(loci[:9]) + "\n")              # 8 loci
    (d / "chunk1.txt").write_text("\n".join(loci[:5]) + "\n")
    (d / "chunk2.txt").write_text("\n".join([loci[0]] + loci[5:9]) + "\n")
    return d


def runner(d, ref, loci, out, info="info.txt"):
    run(["Rscript", BIN / "run_lava.R", "--info", d / info, "--ref", ref, "--loci", d / loci,
         "--pairs", d / "pairs.tsv", "--out", d / out, "--seed", SEED, "--univ-thresh", "0.05"])
    return pd.read_csv(d / f"{out}.lava.tsv", sep="\t", dtype={"locus": str})


def test_runner_matches_native_lava_exactly(lava_inputs, vignette):
    d = lava_inputs
    ours = runner(d, vignette / "g1000_test", "loci.txt", "ours")
    run(["Rscript", ORACLE, d / "info.txt", vignette / "g1000_test", d / "loci.txt",
         d / "native.tsv", SEED, 0.05, *PAIR])
    native = pd.read_csv(d / "native.tsv", sep="\t", dtype={"locus": str})

    m = ours.merge(native, on="locus", suffixes=("", "_native"))
    assert len(m) == len(native) > 0
    tested = m["status"] == "ok"
    assert tested.sum() >= 3, "need some loci that pass the univariate filter"
    for col in ["h2_1", "h2_2", "rg", "rg_ci_low", "rg_ci_high", "p"]:
        np.testing.assert_array_equal(m[col].to_numpy(), m[f"{col}_native"].to_numpy(),
                                      err_msg=f"column {col} differs from native LAVA")


def test_chunking_does_not_change_results(lava_inputs, vignette):
    d = lava_inputs
    whole = runner(d, vignette / "g1000_test", "loci.txt", "whole")
    parts = pd.concat([runner(d, vignette / "g1000_test", "chunk1.txt", "c1"),
                       runner(d, vignette / "g1000_test", "chunk2.txt", "c2")])
    key = ["trait1", "trait2", "locus"]
    pd.testing.assert_frame_equal(whole.sort_values(key).reset_index(drop=True),
                                  parts.sort_values(key).reset_index(drop=True))


def test_harmonized_input_matches_lava_reading_raw_files(lava_inputs, vignette):
    """LAVA's own header handling on the raw files gives the same estimates."""
    d = lava_inputs
    write_info(d / "info_raw.txt", {t: vignette / f"{t}.sumstats.txt" for t in PAIR})
    ours = runner(d, vignette / "g1000_test", "loci.txt", "ours_h")
    raw = runner(d, vignette / "g1000_test", "loci.txt", "raw", info="info_raw.txt")
    m = ours.merge(raw, on=["trait1", "trait2", "locus"], suffixes=("", "_raw"))
    assert (m["status"] == m["status_raw"]).all()
    assert (m["n_snps"] == m["n_snps_raw"]).all()
    # neuro has an explicit Z column: identical input, identical local h2
    np.testing.assert_array_equal(m["h2_1"], m["h2_1_raw"])
    # depression: only the OR == 1 sign convention differs (see module docstring)
    np.testing.assert_allclose(m["h2_2"], m["h2_2_raw"], rtol=1e-3)
    np.testing.assert_allclose(m["rg"], m["rg_raw"], atol=1e-4, equal_nan=True)
