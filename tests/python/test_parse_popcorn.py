"""bin/parse_popcorn.py on real Popcorn output (tests/python/fixtures/popcorn),
made with the gencorr-popcorn image from Popcorn's own EUR/EAS test data after
gencorr harmonization and export."""

import json

import pandas as pd
import pytest
from scipy import stats

from conftest import BIN, FIXTURES, run

import parse_popcorn as pp

POP = FIXTURES / "popcorn"


def pair(stem):
    return pp.parse_pair(POP / f"{stem}.popcorn", POP / f"{stem}.popcorn.o.log",
                         POP / f"{stem}.popcorn.e.log")


def test_values_are_taken_verbatim():
    r = pair("eur_t__eas_t")
    assert (r["trait1"], r["trait2"]) == ("eur_t", "eas_t")
    assert r["rg"] == pytest.approx(0.5195109138379685, abs=0)
    assert r["rg_se"] == pytest.approx(0.15396054562902098, abs=0)
    assert r["h2_1"] == pytest.approx(0.1965337777517921, abs=0)
    assert r["h2_2"] == pytest.approx(0.07685208160292506, abs=0)
    assert r["n_snps"] == 43897 and r["status"] == "ok"


def test_p_column_tests_rg_not_zero_and_popcorn_p_is_kept_separately():
    r = pair("eur_t__eas_t")
    extra = json.loads(r["extra"])
    # standard columns: rg != 0
    assert r["z"] == pytest.approx(0.5195109138379685 / 0.15396054562902098)
    assert r["p"] == pytest.approx(2 * stats.norm.sf(0.5195109138379685 / 0.15396054562902098))
    # Popcorn's own test is rg < 1: Z = (1 - rg) / SE
    assert extra["z_pg_lt_1"] == pytest.approx((1 - 0.5195109138379685) / 0.15396054562902098)
    assert extra["p_pg_lt_1"] == pytest.approx(0.0018032462387601589)
    assert r["p"] != pytest.approx(extra["p_pg_lt_1"])


def test_extra_describes_the_model():
    extra = json.loads(pair("eur_t__eas_t")["extra"])
    assert extra["correlation"] == "genetic_impact" and extra["fit"] == "regression"
    assert extra["popcorn_version"] == "1.1"
    assert extra["h2_1_liab"] is None and extra["error"] is None


def test_liability_scale_columns():
    r = pair("eur_l__eas_l")
    extra = json.loads(r["extra"])
    assert r["rg"] == pytest.approx(0.5195109138379685)          # rg is scale-free
    assert extra["h2_1_liab"] == pytest.approx(0.2067455766236097)
    assert extra["h2_2_liab"] == pytest.approx(0.10458910123851392)


def test_failed_run_is_reported_not_dropped():
    r = pair("empty_t__eas_t")
    assert r["status"] == "popcorn_error" and pd.isna(r["rg"])
    assert "IndexError" in json.loads(r["extra"])["error"]


def test_cli_writes_standard_table(tmp_path):
    run(["python3", BIN / "parse_popcorn.py", "--results", *sorted(POP.glob("*.popcorn")),
         "--logs", *sorted(POP.glob("*.log")), "--out", tmp_path / "gencorr"])
    df = pd.read_csv(tmp_path / "gencorr.popcorn.tsv", sep="\t", quoting=3)
    assert sorted(df.status) == ["ok", "ok", "popcorn_error"]
    assert set(df.method) == {"popcorn"}
