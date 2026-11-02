"""Unit tests for parse_ldsc.py, merge_results.py and compare_results.py."""

import json

import numpy as np
import pandas as pd
import pytest

from conftest import BIN, FIXTURES, run

import parse_ldsc
from gencorr_common import RESULT_COLUMNS, read_tsv, validate_results, write_tsv

LOGS = FIXTURES / "ldsc"


# ---------------------------------------------------------------- parse_ldsc --
def test_h2_values_are_taken_verbatim_from_logs():
    h2 = parse_ldsc.collect_h2([LOGS / "depression.h2.log", LOGS / "depression.h2_liab.log"])
    row = h2.iloc[0]
    # numbers as printed in the fixture logs
    assert (row.h2_obs, row.h2_obs_se) == (0.0023, 0.0003)
    assert (row.h2_liab, row.h2_liab_se) == (0.0031, 0.0004)
    assert (row.intercept, row.intercept_se) == (1.1224, 0.0403)
    # samp_prev = n_cases / (n_cases + n_controls) = 170756 / 500199
    assert (row.samp_prev, row.pop_prev) == pytest.approx((170756 / 500199, 0.15), abs=1e-6)


def test_liability_conversion_matches_ldsc():
    """LDSC's liability h2 equals the Lee et al. 2011 transform of its observed h2
    (up to the 4-decimal rounding of the log), i.e. the two runs are consistent."""
    h2 = parse_ldsc.collect_h2([LOGS / "depression.h2.log", LOGS / "depression.h2_liab.log"]).iloc[0]
    K, P = h2.pop_prev, h2.samp_prev
    t = __import__("scipy.stats", fromlist=["norm"]).norm.isf(K)
    phi = __import__("scipy.stats", fromlist=["norm"]).norm.pdf(t)
    liab = h2.h2_obs * K**2 * (1 - K) ** 2 / (P * (1 - P) * phi**2)
    assert liab == pytest.approx(h2.h2_liab, abs=1e-4)


def test_rg_log_with_several_pairs():
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__b1.rg.log")
    assert list(rg.trait2) == ["neuro", "depression", "rheuma"]
    assert rg.set_index("trait2").loc["neuro", "rg"] == 0.1413
    assert rg.set_index("trait2").loc["neuro", "n_snps"] == 97121
    assert rg.error.isna().all()


def test_rg_error_is_reported_not_dropped():
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__err.rg.log")
    assert len(rg) == 1
    assert "Incompatible alleles" in rg.error.iloc[0]
    h2 = parse_ldsc.collect_h2([LOGS / "bmi.h2.log", LOGS / "depression.h2.log"])
    res = parse_ldsc.build_results(rg, h2)
    assert res.status.iloc[0] == "ldsc_error"
    assert np.isnan(res.rg.iloc[0])


def test_sample_overlap_matrix_follows_lava_vignette():
    h2 = parse_ldsc.collect_h2([LOGS / f"{t}.h2.log" for t in ["bmi", "neuro", "depression", "rheuma"]])
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__b1.rg.log")
    mat = parse_ldsc.sample_overlap_matrix(rg, h2)
    assert (np.diag(mat.values) == 1).all()
    assert np.allclose(mat.values, mat.values.T)
    hi = h2.set_index("trait").intercept
    expected = 0.0158 / np.sqrt(hi["bmi"] * hi["neuro"])      # gcov_int / sqrt(int1 * int2)
    assert mat.loc["bmi", "neuro"] == pytest.approx(round(expected, 5))
    assert mat.loc["neuro", "rheuma"] == 0.0                   # pair not run -> no correction


def test_parse_ldsc_cli(tmp_path):
    logs = [LOGS / f"{t}.h2.log" for t in ["bmi", "neuro", "depression", "rheuma"]]
    run(["python3", BIN / "parse_ldsc.py", "--h2-logs", *logs, LOGS / "depression.h2_liab.log",
         "--rg-logs", LOGS / "bmi__b1.rg.log", "--out", tmp_path / "x"])
    res = validate_results(read_tsv(tmp_path / "x.ldsc.tsv"))
    assert len(res) == 3 and set(res.status) == {"ok"}
    assert json.loads(res.extra.iloc[0])["gcov_int"] is not None
    # overlap file is readable the way LAVA reads it (header lacks the row-name column)
    lines = (tmp_path / "x.sample_overlap.txt").read_text().splitlines()
    assert len(lines[0].split()) + 1 == len(lines[1].split())


# --------------------------------------------------------------------- merge --
def _results(**over):
    row = {c: np.nan for c in RESULT_COLUMNS}
    row.update(method="ldsc", scope="global", trait1="a", trait2="b", rg=0.5, status="ok", extra="{}")
    row.update(over)
    return pd.DataFrame([row])


def test_merge_rejects_duplicates(tmp_path):
    write_tsv(_results(), tmp_path / "1.tsv")
    write_tsv(_results(), tmp_path / "2.tsv")
    with pytest.raises(AssertionError, match="duplicated"):
        run(["python3", BIN / "merge_results.py", "--out", tmp_path / "m.tsv",
             tmp_path / "1.tsv", tmp_path / "2.tsv"])


def test_merge_rejects_wrong_schema(tmp_path):
    bad = _results().rename(columns={"rg": "rho"})
    bad.to_csv(tmp_path / "bad.tsv", sep="\t", index=False)
    with pytest.raises(AssertionError, match="missing result columns"):
        run(["python3", BIN / "merge_results.py", "--out", tmp_path / "m.tsv", tmp_path / "bad.tsv"])


# ------------------------------------------------------------------- compare --
def test_compare_pass_fail_and_order_insensitive(tmp_path):
    write_tsv(_results(trait1="a", trait2="b", rg=0.5), tmp_path / "test.tsv")
    write_tsv(_results(trait1="b", trait2="a", rg=0.50005), tmp_path / "ref.tsv")
    run(["python3", BIN / "compare_results.py", "--test", tmp_path / "test.tsv",
         "--reference", tmp_path / "ref.tsv", "--method", "ldsc"])
    write_tsv(_results(trait1="b", trait2="a", rg=0.6), tmp_path / "ref2.tsv")
    with pytest.raises(AssertionError, match="FAIL"):
        run(["python3", BIN / "compare_results.py", "--test", tmp_path / "test.tsv",
             "--reference", tmp_path / "ref2.tsv", "--method", "ldsc"])


def test_compare_reads_bigagwas_ldsc_csv(tmp_path):
    write_tsv(pd.concat([_results(trait1="t1", trait2="x", rg=0.25, rg_se=0.05, z=5.0, p=5.7e-7),
                         _results(trait1="t1", trait2="y", rg=np.nan, status="rg_not_estimable")]),
              tmp_path / "test.tsv")
    biga = pd.DataFrame({
        "Trait1": ["MyTrait", "MyTrait"], "Trait2": ["x", "y"], "Trait2_name": ["X", "Y"],
        "group": ["g", "g"], "rg": ["0.2500", "9999"], "se": ["0.0500", "9999"],
        "z": ["5.0000", "9999"], "p_value": ["5.7000e-07", "9999"]})
    biga.to_csv(tmp_path / "biga.csv", index=False)
    pd.DataFrame({"reference_name": ["MyTrait"], "gencorr_id": ["t1"]}).to_csv(tmp_path / "map.csv", index=False)
    run(["python3", BIN / "compare_results.py", "--test", tmp_path / "test.tsv", "--method", "ldsc",
         "--reference", tmp_path / "biga.csv", "--reference-format", "bigagwas_ldsc",
         "--trait-map", tmp_path / "map.csv", "--report", tmp_path / "rep"])
    assert "matched pairs      : 2" in (tmp_path / "rep.txt").read_text()
