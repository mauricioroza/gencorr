"""bin/parse_sumher.py on real LDAK 6.1 output (tests/python/fixtures/sumher)."""

import json

import pandas as pd
import pytest
from scipy import stats

from conftest import BIN, FIXTURES, run

import parse_sumher as ps

SUMHER = FIXTURES / "sumher"


def test_values_are_taken_verbatim_from_cors():
    r = ps.parse_pair(SUMHER / "bmi__depression.log", SUMHER / "bmi__depression.cors")
    assert (r["trait1"], r["trait2"]) == ("bmi", "depression")
    # numbers as printed in the fixture .cors
    assert (r["rg"], r["rg_se"]) == (-0.093256, 0.057574)
    assert (r["h2_1"], r["h2_1_se"]) == (0.004497, 0.000385)
    assert (r["h2_2"], r["h2_2_se"]) == (0.001186, 0.000094)
    assert r["status"] == "ok"


def test_z_p_ci_are_normal_approximation():
    r = ps.parse_pair(SUMHER / "bmi__depression.log", SUMHER / "bmi__depression.cors")
    assert r["z"] == pytest.approx(-0.093256 / 0.057574)
    assert r["p"] == pytest.approx(2 * stats.norm.sf(0.093256 / 0.057574))
    assert (r["rg_ci_low"], r["rg_ci_high"]) == pytest.approx((-0.093256 - 1.96 * 0.057574,
                                                               -0.093256 + 1.96 * 0.057574))


def test_extra_has_coheritability_overlap_and_snp_counts():
    extra = json.loads(ps.parse_pair(SUMHER / "bmi__depression.log", SUMHER / "bmi__depression.cors")["extra"])
    assert (extra["coher"], extra["coher_se"]) == (-0.000215, 0.000131)
    assert (extra["overlap"], extra["overlap_se"]) == (0.223522, 0.045501)
    assert (extra["n_snps_1"], extra["n_snps_2"]) == (39310, 28538)
    assert "scaling_1" not in extra          # --genomic-control NO: no scaling rows
    assert extra["error"] is None


def test_scaling_rows_are_kept_when_present(tmp_path):
    cors = tmp_path / "x.cors"
    cors.write_text("Component Value SE\nHer1_All 0.1 0.01\nHer2_All 0.2 0.02\nCoher_All 0.05 0.01\n"
                    "Cor_All 0.35 0.05\nScaling1 1.1 0.02\nScaling2 0.9 0.03\nOverlap 0.1 0.02\n")
    extra = json.loads(ps.parse_pair(SUMHER / "bmi__depression.log", cors)["extra"])
    assert (extra["scaling_1"], extra["scaling_2"]) == (1.1, 0.9)


def test_failed_ldak_run_is_reported_not_dropped():
    r = ps.parse_pair(SUMHER / "empty__depression.log", SUMHER / "empty__depression.cors")
    assert r["status"] == "sumher_error"
    assert pd.isna(r["rg"]) and pd.isna(r["p"])
    assert "contains summary statistics for none" in json.loads(r["extra"])["error"]


def test_cli_writes_standard_table(tmp_path):
    run(["python3", BIN / "parse_sumher.py",
         "--logs", SUMHER / "bmi__depression.log", SUMHER / "empty__depression.log",
         "--cors", SUMHER / "bmi__depression.cors", SUMHER / "empty__depression.cors",
         "--out", tmp_path / "gencorr"])
    df = pd.read_csv(tmp_path / "gencorr.sumher.tsv", sep="\t", quoting=3)
    assert list(df.status) == ["ok", "sumher_error"]
    assert set(df.method) == {"sumher"} and set(df.scope) == {"global"}
