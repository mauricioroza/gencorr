"""Unit tests for bin/harmonize_sumstats.py and bin/export_sumstats.py."""

import gzip
import hashlib
import json

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from conftest import BIN, run

import harmonize_sumstats as hs
import export_sumstats as ex


def raw(**cols):
    """Build a raw (all-string) table like read_table() returns."""
    return pd.DataFrame({k: [None if v is None else str(v) for v in vals] for k, vals in cols.items()})


def harmonize(df, overrides=None, **kw):
    mapping = hs.map_columns(list(df.columns), overrides or {})
    return hs.harmonize(df, mapping, **kw)


# ----------------------------------------------------------------- mapping --
def test_column_synonyms_are_case_insensitive():
    m = hs.map_columns(["MarkerName", "Allele1", "Allele2", "Zscore", "P-value", "TotalN"], {})
    assert m == {"SNP": "MarkerName", "A1": "Allele1", "A2": "Allele2", "Z": "Zscore", "P": "P-value"}


def test_stat_is_not_mapped_automatically():
    # STAT may be a t, Z or chi2 depending on the software: must be explicit
    m = hs.map_columns(["SNP", "A1", "A2", "STAT"], {})
    assert "Z" not in m
    assert hs.map_columns(["SNP", "A1", "A2", "STAT"], {"Z": "STAT"})["Z"] == "STAT"


def test_override_unknown_column_fails():
    with pytest.raises(SystemExit):
        hs.map_columns(["SNP", "A1", "A2"], {"Z": "nope"})


# ------------------------------------------------------------------ Z rules --
def test_z_priority_explicit_z_first():
    df = raw(SNP=["rs1"], A1=["a"], A2=["g"], Z=["2.5"], BETA=["-1"], SE=["1"], P=["0.5"], N=["100"])
    out, qc = harmonize(df)
    assert out.loc[0, "Z"] == 2.5
    assert qc["z_source_counts"] == {"z": 1}


def test_z_from_p_and_beta_sign_matches_ldsc_convention():
    df = raw(SNP=["rs1", "rs2"], A1=["A", "A"], A2=["G", "G"], BETA=["-0.2", "0.1"],
             SE=["0.1", "0.1"], P=["0.01", "0.2"], N=["1000", "1000"])
    out, qc = harmonize(df)
    expected = np.array([-stats.norm.isf(0.005), stats.norm.isf(0.1)])
    assert np.allclose(out["Z"], expected, rtol=0, atol=1e-12)
    assert qc["z_source_counts"] == {"p": 2}


def test_z_falls_back_to_beta_se_when_p_unusable():
    df = raw(SNP=["rs1"], A1=["A"], A2=["G"], BETA=["0.5"], SE=["0.1"], P=["0"], N=["50"])
    out, qc = harmonize(df)
    assert out.loc[0, "Z"] == pytest.approx(5.0)
    assert qc["z_source_counts"] == {"beta_se": 1}
    # P is recomputed from Z, not left at 0
    assert out.loc[0, "P"] == pytest.approx(2 * stats.norm.sf(5.0))


def test_or_only_gives_log_or_beta_and_sign():
    df = raw(SNP=["rs1", "rs2"], A1=["A", "A"], A2=["G", "G"], OR=["0.8", "1.25"],
             P=["0.05", "0.05"], N=["10", "10"])
    out, qc = harmonize(df)
    assert np.allclose(out["BETA"], np.log([0.8, 1.25]))
    assert (np.sign(out["Z"]) == [-1, 1]).all()
    assert qc["beta_from_or"] is True


def test_z_priority_can_be_changed():
    df = raw(SNP=["rs1"], A1=["A"], A2=["G"], BETA=["0.5"], SE=["0.1"], P=["0.5"], N=["50"])
    out, _ = harmonize(df, z_priority=("beta_se", "p"))
    assert out.loc[0, "Z"] == pytest.approx(5.0)


# ------------------------------------------------------------------ filters --
def test_filters_are_counted():
    df = raw(
        SNP=["rs1", "rs2", "rs3", "rs4", "rs4", "rs5", "RS6"],
        A1=["A", "AT", "C", "A", "A", "T", "c"],
        A2=["G", "A", "C", "G", "G", "A", "g"],
        Z=["1", "1", "1", "1", "1", "1", "1"],
        N=["10", "10", "10", "10", "10", "10", "10"],
    )
    out, qc = harmonize(df, drop_ambiguous=True)
    assert list(out["SNP"]) == ["rs1"]            # rs6 is C/G (ambiguous) after upper-casing
    assert qc["removed"]["non_snv_or_missing_allele"] == 1
    assert qc["removed"]["a1_equals_a2"] == 1
    assert qc["removed"]["strand_ambiguous"] == 2     # rs5 A/T, rs6 C/G
    assert qc["removed"]["duplicated_snp_all_copies"] == 2
    assert qc["rows_in"] == 7 and qc["rows_out"] == 1


def test_rs_ids_lower_cased_and_alleles_upper_cased():
    df = raw(SNP=["RS123", "chr1:100"], A1=["a", "c"], A2=["g", "t"], Z=["1", "2"], N=["5", "5"])
    out, _ = harmonize(df)
    assert list(out["SNP"]) == ["rs123", "chr1:100"]
    assert list(out["A1"]) == ["A", "C"]


def test_n_fallbacks():
    df = raw(SNP=["rs1"], A1=["A"], A2=["G"], Z=["1"])
    out, qc = harmonize(df, n_cases=100, n_controls=300)
    assert out.loc[0, "N"] == 400 and qc["n_source"] == "samplesheet"
    with pytest.raises(SystemExit):
        harmonize(df)


def test_effective_n_column_is_flagged():
    df = raw(SNP=["rs1"], A1=["A"], A2=["G"], Z=["1"], Neff=["123"])
    _, qc = harmonize(df)
    assert "warning_effective_n" in qc


# ------------------------------------------------------------- CLI + export --
def test_cli_output_is_byte_identical_across_runs(tmp_path):
    src = tmp_path / "in.tsv"
    pd.DataFrame({"SNP": ["rs1", "rs2"], "A1": ["A", "C"], "A2": ["G", "T"],
                  "BETA": [0.1, -0.2], "SE": [0.05, 0.05], "P": [0.04, 1e-4],
                  "N": [1000, 1000]}).to_csv(src, sep="\t", index=False)
    digests = []
    for i in range(2):                     # same output name, as in the pipeline
        (tmp_path / f"run{i}").mkdir()
        out = tmp_path / f"run{i}" / "t"
        run(["python3", BIN / "harmonize_sumstats.py", "--input", src, "--id", "t", "--out", out])
        digests.append(hashlib.md5((tmp_path / f"run{i}" / "t.harmonized.tsv.gz").read_bytes()).hexdigest())
    assert digests[0] == digests[1]
    qc = json.loads((tmp_path / "run0" / "t.qc.json").read_text())
    assert qc["rows_out"] == 2


def test_export_formats_and_optional_columns(tmp_path):
    h = pd.DataFrame({"SNP": ["rs1", "rs2"], "CHR": ["1", "1"], "BP": [1, 2], "A1": ["A", "C"],
                      "A2": ["G", "T"], "Z": [1.0, -2.0], "BETA": [np.nan] * 2, "SE": [np.nan] * 2,
                      "P": [0.3, 0.04], "N": [10.0, 10.0], "EAF": [0.2, np.nan],
                      "INFO": [0.95, 0.99]})
    ldsc = ex.export(h, "ldsc")
    # INFO complete -> kept; EAF half missing -> dropped (LDSC would drop those rows)
    assert list(ldsc.columns) == ["SNP", "A1", "A2", "Z", "P", "N", "INFO"]
    lava = ex.export(h, "lava")
    assert list(lava.columns) == ["SNP", "A1", "A2", "Z", "N"]

    src = tmp_path / "t.harmonized.tsv.gz"
    h.to_csv(src, sep="\t", index=False, na_rep="NA", compression="gzip")
    run(["python3", BIN / "export_sumstats.py", "--input", src, "--format", "ldsc",
         "--out", tmp_path / "t.ldsc.tsv"])
    with gzip.open(tmp_path / "t.ldsc.tsv.gz", "rt") as fh:
        assert fh.readline().split() == ["SNP", "A1", "A2", "Z", "P", "N", "INFO"]


@pytest.mark.parametrize("trait,expected_source", [
    ("bmi", "p"), ("neuro", "z"), ("depression", "p"), ("rheuma", "p")])
def test_real_vignette_layouts(vignette, tmp_path, trait, expected_source):
    """Every LAVA vignette layout is read without manual column mapping."""
    run(["python3", BIN / "harmonize_sumstats.py", "--input", vignette / f"{trait}.sumstats.txt",
         "--id", trait, "--out", tmp_path / trait])
    qc = json.loads((tmp_path / f"{trait}.qc.json").read_text())
    assert qc["rows_out"] > 50_000
    assert set(qc["z_source_counts"]) == {expected_source}
