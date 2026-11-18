"""Allele orientation and column detection for the output layouts of common
GWAS software (bin/harmonize_sumstats.py: COLUMN_SYNONYMS, detect_layout).

Each test writes a few rows in the software's own header layout and reads them
through load_input(), as the pipeline does."""

import numpy as np
import pytest
from scipy import stats

import harmonize_sumstats as hs


def harmonize_file(tmp_path, lines, overrides=None, **kw):
    src = tmp_path / "in.txt"
    src.write_text("\n".join("\t".join(r.split()) for r in lines) + "\n")
    raw, mapping = hs.load_input(src, overrides or {})
    out, qc = hs.harmonize(raw, mapping, **kw)
    return out.set_index("SNP"), qc


SAIGE = [
    "CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 imputationInfo N BETA SE Tstat var p.value",
    "1 100 rs1 A G 500 0.25 0.95 1000 0.2 0.08 31.2 156 0.0124",
    "1 200 rs2 C T 300 0.15 0.90 1000 -0.1 0.05 -40 400 0.0455",
]


def test_saige_effect_allele_is_allele2(tmp_path):
    out, qc = harmonize_file(tmp_path, SAIGE)
    assert qc["layout"] == "saige"
    assert (out.loc["rs1", "A1"], out.loc["rs1", "A2"]) == ("G", "A")
    assert (out.loc["rs2", "A1"], out.loc["rs2", "A2"]) == ("T", "C")
    assert out.loc["rs1", "Z"] > 0 and out.loc["rs2", "Z"] < 0      # sign follows BETA of Allele2
    assert out.loc["rs1", "EAF"] == 0.25                              # AF_Allele2, same orientation
    assert out.loc["rs1", "INFO"] == 0.95
    assert "warning_allele_orientation" not in qc


def test_colmap_alleles_switch_layout_off(tmp_path):
    out, qc = harmonize_file(tmp_path, SAIGE, overrides={"A1": "Allele1", "A2": "Allele2"})
    assert "layout" not in qc
    assert out.loc["rs1", "A1"] == "A"


PLINK2_HEADER = "#CHROM POS ID REF ALT A1 A1_FREQ OBS_CT BETA SE T_STAT P"
PLINK2_ROWS = [
    "1 100 rs1 A G G 0.3 1000 0.2 0.05 4 6.3e-05",    # ALT tested
    "1 200 rs2 C T C 0.6 1000 -0.1 0.05 -2 0.0455",   # REF tested
    "1 300 rs3 G A A 0.2 1000 0.05 0.05 1 0.317",
]


def test_plink2_a1_equal_to_ref_keeps_variant(tmp_path):
    out, qc = harmonize_file(tmp_path, [PLINK2_HEADER] + PLINK2_ROWS)
    assert qc["layout"] == "plink2" and qc["a2_from_alt"] == 1
    assert qc["removed"]["a1_equals_a2"] == 0 and qc["rows_out"] == 3
    assert [tuple(out.loc[s, ["A1", "A2"]]) for s in ("rs1", "rs2", "rs3")] == [
        ("G", "A"), ("C", "T"), ("A", "G")]
    assert out.loc["rs2", "Z"] < 0
    assert out.loc["rs2", "EAF"] == 0.6                               # A1_FREQ


def test_plink2_omitted_column_is_a2(tmp_path):
    header = PLINK2_HEADER.replace("A1 A1_FREQ", "A1 OMITTED A1_FREQ")
    rows = [r.split() for r in PLINK2_ROWS]
    omitted = ["A", "T", "G"]
    lines = [header] + [" ".join(r[:6] + [o] + r[6:]) for r, o in zip(rows, omitted)]
    out, qc = harmonize_file(tmp_path, lines)
    assert qc["layout"] == "plink2" and "a2_from_alt" not in qc
    assert list(out["A2"]) == omitted


def test_plink2_log10_p(tmp_path):
    lines = [PLINK2_HEADER.rsplit(" ", 1)[0] + " LOG10_P"] + [
        " ".join(r.split()[:-1] + [str(-np.log10(float(r.split()[-1])))]) for r in PLINK2_ROWS]
    out, qc = harmonize_file(tmp_path, lines)
    assert qc["p_from_log10p"] is True
    assert out.loc["rs2", "P"] == pytest.approx(0.0455)


def test_regenie(tmp_path):
    lines = [
        "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P",
        f"1 12345 rs1 A G 0.3 0.97 5000 ADD 0.1 0.05 4 {-np.log10(0.05)}",
    ]
    out, qc = harmonize_file(tmp_path, lines)
    assert (out.loc["rs1", "A1"], out.loc["rs1", "A2"]) == ("G", "A")
    assert out.loc["rs1", "BP"] == 12345
    assert qc["z_source_counts"] == {"p": 1}
    assert out.loc["rs1", "Z"] == pytest.approx(stats.norm.isf(0.025))
    assert "warning_allele_orientation" not in qc


def test_bolt_lmm(tmp_path):
    lines = [
        "SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ F_MISS BETA SE P_BOLT_LMM_INF P_BOLT_LMM",
        "rs1 1 12345 0.5 G A 0.3 0.01 -0.1 0.05 0.04 0.05",
    ]
    out, qc = harmonize_file(tmp_path, lines, n=1000)
    assert out.loc["rs1", "BP"] == 12345                               # BP, not genetic GENPOS
    assert (out.loc["rs1", "A1"], out.loc["rs1", "A2"]) == ("G", "A")
    assert qc["column_mapping"]["P"] == "P_BOLT_LMM"
    assert out.loc["rs1", "Z"] == pytest.approx(-stats.norm.isf(0.025))


def test_metal_allele1_is_effect_without_warning(tmp_path):
    lines = ["MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value TotalSampleSize",
             "rs1 a g 0.3 0.1 0.05 0.05 1000"]
    out, qc = harmonize_file(tmp_path, lines, n=1000)
    assert out.loc["rs1", "A1"] == "A" and out.loc["rs1", "Z"] > 0
    assert "warning_allele_orientation" not in qc and "layout" not in qc


def test_unknown_allele1_allele2_layout_is_flagged(tmp_path):
    lines = ["SNP Allele1 Allele2 BETA SE P N", "rs1 A G 0.1 0.05 0.05 1000"]
    _, qc = harmonize_file(tmp_path, lines)
    assert "warning_allele_orientation" in qc


def test_r2_is_not_taken_as_info():
    assert "INFO" not in hs.map_columns(["SNP", "A1", "A2", "Z", "R2"], {})
