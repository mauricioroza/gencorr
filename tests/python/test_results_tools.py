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


def test_h2_table_does_not_depend_on_log_order():
    """Nextflow stages and globs files in no guaranteed order; reading the
    liability log first used to blank samp_prev/pop_prev."""
    logs = [LOGS / "depression.h2.log", LOGS / "depression.h2_liab.log"]
    fwd = parse_ldsc.collect_h2(logs)
    rev = parse_ldsc.collect_h2(logs[::-1])
    pd.testing.assert_frame_equal(fwd, rev)
    assert rev.iloc[0].samp_prev == pytest.approx(170756 / 500199, abs=1e-6)


def test_h2_status_reports_failed_and_empty_logs(tmp_path, capsys):
    """A trait whose LDSC h2 run failed keeps its row, with a status and a warning."""
    good = (LOGS / "bmi.h2.log").read_text()
    (tmp_path / "bmi.h2.log").write_text(good)
    (tmp_path / "bad.h2.log").write_text(
        "* Version 1.0.1\nTraceback (most recent call last):\n  File \"ldsc.py\"\n"
        "ValueError: After merging with reference panel LD, 0 SNPs remain.\n")
    (tmp_path / "empty.h2.log").write_text("* Version 1.0.1\nAnalysis finished\n")
    h2 = parse_ldsc.collect_h2(sorted(tmp_path.glob("*.h2.log"))).set_index("trait")
    assert dict(h2.status) == {"bad": "h2_error", "bmi": "ok", "empty": "h2_not_estimable"}
    assert np.isnan(h2.loc["bad", "h2_obs"])
    err = capsys.readouterr().err
    assert "LDSC h2 failed for 'bad': ValueError: After merging" in err
    assert "no observed-scale h2 in the LDSC log of 'empty'" in err


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


def test_sample_overlap_is_zero_for_trait_without_intercept():
    """A trait whose h2 intercept is missing gets no overlap correction at all:
    its gcov_int values are not scaled with a made-up diagonal."""
    h2 = parse_ldsc.collect_h2([LOGS / f"{t}.h2.log" for t in ["bmi", "neuro", "depression", "rheuma"]])
    h2.loc[h2.trait == "neuro", "intercept"] = np.nan
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__b1.rg.log")
    mat = parse_ldsc.sample_overlap_matrix(rg, h2)
    assert mat.loc["bmi", "neuro"] == 0.0 and mat.loc["neuro", "bmi"] == 0.0
    assert mat.loc["bmi", "depression"] != 0.0                # other pairs unchanged
    assert (np.diag(mat.values) == 1).all()


def test_overlap_problems_reports_invalid_matrices():
    h2 = parse_ldsc.collect_h2([LOGS / f"{t}.h2.log" for t in ["bmi", "neuro", "depression", "rheuma"]])
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__b1.rg.log")
    assert parse_ldsc.overlap_problems(parse_ldsc.sample_overlap_matrix(rg, h2)) == []
    assert parse_ldsc.overlap_problems(pd.DataFrame()) == []

    t = ["a", "b", "c"]
    too_big = pd.DataFrame([[1, 1.2, 0], [1.2, 1, 0], [0, 0, 1]], index=t, columns=t)
    problems = parse_ldsc.overlap_problems(too_big)
    assert any("a x b" in p for p in problems) and any("eigenvalue" in p for p in problems)
    # every pair valid (|r| < 1) but the three together are not: only joint mode is affected
    joint_only = pd.DataFrame([[1, 0.9, -0.9], [0.9, 1, 0.9], [-0.9, 0.9, 1]], index=t, columns=t)
    problems = parse_ldsc.overlap_problems(joint_only)
    assert len(problems) == 1 and "joint" in problems[0]


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


def test_ldsc_extra_has_diagnostics_and_overlap_flag():
    h2 = parse_ldsc.collect_h2([LOGS / f"{t}.h2.log" for t in ["bmi", "neuro", "depression", "rheuma"]])
    rg = parse_ldsc.parse_rg_log(LOGS / "bmi__b1.rg.log")
    res = parse_ldsc.build_results(rg, h2).set_index("trait2")
    hi = h2.set_index("trait")
    x = json.loads(res.loc["neuro", "extra"])
    assert x["mean_chi2_1"] == hi.loc["bmi", "mean_chi2"] and x["mean_chi2_2"] == hi.loc["neuro", "mean_chi2"]
    assert x["intercept_1"] == hi.loc["bmi", "intercept"]
    assert x["overlap_suspected"] is (abs(x["gcov_int"]) > parse_ldsc.OVERLAP_GCOV_INT)
    assert set(res["confidence"].dropna()) <= {"trusted", "low", "not_interpretable"}


# ---------------------------------------------------------------- confidence --
@pytest.mark.parametrize("over, expected", [
    ({}, "trusted"),                                                    # h2 z = 10
    ({"h2_2": 0.03}, "low"),                                            # h2 z = 3
    ({"h2_1": 0.02}, "not_interpretable"),                              # h2 z = 2
    ({"rg": 1.05}, "not_interpretable"),                                # out of bounds
    ({"rg_se": 0.3}, "not_interpretable"),                              # PGC Popcorn SE rule
    ({"extra": '{"mean_chi2_1": 1.01, "mean_chi2_2": 1.3}'}, "low"),    # LDSC FAQ mean chi2
    ({"extra": '{"mean_chi2_1": null}'}, "trusted"),                    # missing chi2 ignored
    ({"n_snps": 18_000}, "trusted"),                                    # absolute count not used
    ({"extra": '{"snp_fraction_1": 0.85, "snp_fraction_2": 0.13}'}, "not_interpretable"),  # exome-imputed
    ({"extra": '{"snp_fraction_1": 0.85, "snp_fraction_2": 0.6}'}, "trusted"),   # evenly sparser array
    ({"extra": '{"snp_fraction_1": null, "snp_fraction_2": 0.9}'}, "trusted"),   # missing value ignored
    ({"extra": '{"block_coverage_1": 0.99, "block_coverage_2": 0.7}'}, "low"),  # regions missing
    ({"extra": '{"block_coverage_1": 0.99, "block_coverage_2": 0.95}'}, "trusted"),
    ({"h2_1_se": np.nan}, None),                                        # cannot be assessed
    ({"status": "rg_not_estimable"}, None),
    ({"scope": "local", "locus": "1"}, None),
])
def test_confidence_tiers(over, expected):
    from gencorr_common import row_confidence
    row = {"scope": "global", "status": "ok", "rg": 0.4, "rg_se": 0.1,
           "h2_1": 0.1, "h2_1_se": 0.01, "h2_2": 0.1, "h2_2_se": 0.01, "extra": "{}"}
    row.update(over)
    assert row_confidence(row) == expected


def test_validate_results_fills_confidence_for_tables_without_it():
    """Method tables written without the column (e.g. run_lava.R) still validate."""
    df = _results(rg_se=0.1, h2_1=0.1, h2_1_se=0.01, h2_2=0.1, h2_2_se=0.01).drop(columns="confidence")
    out = validate_results(df)
    assert list(out.columns) == RESULT_COLUMNS and out["confidence"].iloc[0] == "trusted"


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


def test_bh_qvalues_match_r_p_adjust():
    import merge_results
    p = np.array([0.01, 0.04, 0.03, 0.5, 0.001, 0.04])
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.5, 0.001, 0.04), "BH")
    expected = [0.03, 0.048, 0.048, 0.5, 0.006, 0.048]
    assert np.allclose(merge_results.bh_qvalues(p), expected)


def test_merge_adjust_p_counts_only_tested_rows(tmp_path):
    lava = dict(method="lava", scope="local", trait1="a", trait2="b", chr="1", start=1, stop=2)
    rows = [
        _results(**lava, locus="1", p=0.01, extra='{"K":2,"univ_p_1":1e-04}'),
        _results(**lava, locus="2", p=0.04, extra='{"K":3}'),
        _results(**lava, locus="3", p=np.nan, status="univ_filtered", extra='{"K":NA}'.replace("NA", "null")),
        _results(**lava, locus="4", p=0.5, status="estimate_out_of_bounds", extra="{}"),
        _results(locus=np.nan, p=0.2),                      # ldsc: its own family of one
    ]
    for i, r in enumerate(rows):
        write_tsv(r, tmp_path / f"{i}.tsv")
    run(["python3", BIN / "merge_results.py", "--adjust-p", "--out", tmp_path / "m.tsv",
         *[tmp_path / f"{i}.tsv" for i in range(len(rows))]])
    m = read_tsv(tmp_path / "m.tsv").set_index(["method", "locus"], drop=False)
    e1 = m.loc[("lava", "1"), "extra"]
    assert e1.startswith('{"K":2,"univ_p_1":1e-04,')                  # R's text kept as is
    x1, x2 = json.loads(e1), json.loads(m.loc[("lava", "2"), "extra"])
    assert (x1["n_tests"], x1["p_bonferroni"], x1["q_fdr"]) == (2, 0.02, 0.02)
    assert (x2["p_bonferroni"], x2["q_fdr"]) == (0.08, 0.04)
    assert "q_fdr" not in json.loads(m.loc[("lava", "3"), "extra"])
    assert json.loads(m.loc[("lava", "4"), "extra"]) == {}
    ldsc = json.loads(m[m["method"] == "ldsc"]["extra"].iloc[0])
    assert ldsc == {"p_bonferroni": 0.2, "q_fdr": 0.2, "n_tests": 1}

    # running the merge again on its own output does not stack the keys
    run(["python3", BIN / "merge_results.py", "--adjust-p", "--out", tmp_path / "m2.tsv", tmp_path / "m.tsv"])
    assert read_tsv(tmp_path / "m2.tsv")["extra"].tolist() == read_tsv(tmp_path / "m.tsv")["extra"].tolist()


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


def test_rg_p_from_python3_ldsc_log_equals_original():
    # CBIIT/ldsc (Python 3, current pandas) prints the summary-table p with 4
    # fixed decimals (1.1824e-34 -> 0.0000); the parser then takes the pair's
    # own "P:" line. Every other value is printed identically by both versions.
    # Same run (neuro vs depression, rheuma) logged by LDSC 1.0.1 and by the port.
    versions = FIXTURES / "ldsc_versions"
    orig = parse_ldsc.parse_rg_log(versions / "neuro__ldsc1.0.1.rg.log").set_index("trait2")
    py3 = parse_ldsc.parse_rg_log(versions / "neuro__ldsc3.0.1.rg.log").set_index("trait2")
    assert py3.loc["depression", "p"] == orig.loc["depression", "p"] == 1.1824e-34
    assert py3.loc["rheuma", "p"] == 0.051           # "P: 0.051" (4 significant digits)
    assert orig.loc["rheuma", "p"] == 0.051006       # LDSC 1.0.1 table: 5.1006e-02
    cols = ["rg", "se", "z", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se", "n_snps"]
    pd.testing.assert_frame_equal(py3[cols], orig[cols])
