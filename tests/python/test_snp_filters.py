"""
Tests for the two optional SNP filters:

  --min_info / --min_maf   QC thresholds in harmonize_sumstats.py (all methods)
  --lava_snp_list          SNP-ID list applied when exporting LAVA input
                           (export_sumstats.py --snp-list)
"""

import gzip

import numpy as np
import pandas as pd
import pytest

from conftest import BIN, run

import export_sumstats as ex
import harmonize_sumstats as hs


def raw(**cols):
    return pd.DataFrame({k: [None if v is None else str(v) for v in vals] for k, vals in cols.items()})


def harmonize(df, **kw):
    return hs.harmonize(df, hs.map_columns(list(df.columns), {}), **kw)


# ------------------------------------------------------------ thresholds --- #
def base(**extra):
    cols = dict(SNP=["rs1", "rs2", "rs3", "rs4"], A1=["A"] * 4, A2=["G"] * 4,
                Z=["1", "2", "3", "4"], N=["100"] * 4)
    cols.update(extra)
    return raw(**cols)


def test_thresholds_off_by_default():
    out, qc = harmonize(base(INFO=["0.1", "0.5", "0.95", None]))
    assert len(out) == 4
    assert qc["thresholds"] == {"min_info": None, "min_maf": None}


def test_min_info_uses_ldsc_rule_and_keeps_missing():
    # LDSC munge removes INFO <= threshold; a missing INFO is kept
    out, qc = harmonize(base(INFO=["0.85", "0.9", "0.95", None]), min_info=0.9)
    assert list(out["SNP"]) == ["rs3", "rs4"]
    assert qc["removed"]["info_le_0.9"] == 2
    assert qc["info_missing_kept"] == 1


def test_min_maf_from_maf_column_or_eaf():
    out, qc = harmonize(base(MAF=["0.005", "0.01", "0.2", "0.3"]), min_maf=0.01)
    assert list(out["SNP"]) == ["rs3", "rs4"] and qc["maf_source"] == "column:MAF"
    # EAF 0.995 is MAF 0.005 (folded); EAF 0.5 stays
    out, qc = harmonize(base(EAF=["0.995", "0.5", "0.02", "0.009"]), min_maf=0.01)
    assert list(out["SNP"]) == ["rs2", "rs3"] and qc["maf_source"] == "EAF:EAF"


def test_threshold_without_column_warns_and_removes_nothing():
    out, qc = harmonize(base(), min_info=0.9, min_maf=0.01)
    assert len(out) == 4
    assert "warning_min_info" in qc and "warning_min_maf" in qc


def test_maf_column_is_not_written_to_harmonized_file():
    out, _ = harmonize(base(MAF=["0.2"] * 4), min_maf=0.01)
    assert "MAF" not in out.columns and "_MAF" not in out.columns


# -------------------------------------------------------------- SNP list --- #
def test_read_snp_list_header_case_and_gzip(tmp_path):
    lst = tmp_path / "w_hm3.snplist.gz"
    with gzip.open(lst, "wt") as fh:
        fh.write("SNP A1 A2\nRS1 A G\nrs3 C T\n\nchr1:500 A C\n")
    assert ex.read_snp_list(str(lst)) == {"rs1", "rs3", "chr1:500"}


def test_snp_list_cli_equals_manual_awk_filter(tmp_path):
    """--snp-list gives the same LAVA input as filtering the raw file by ID first
    (a manual `awk` step)."""
    snps = [f"rs{i}" for i in range(1, 11)]
    raw_df = pd.DataFrame({"SNP": snps, "A1": ["A"] * 10, "A2": ["G"] * 10,
                           "BETA": np.linspace(-0.1, 0.1, 10), "SE": [0.05] * 10,
                           "P": np.linspace(0.01, 0.9, 10), "N": [1000] * 10})
    keep = {"rs2", "rs3", "rs7", "rs99"}
    (tmp_path / "list.txt").write_text("SNP A1 A2\n" + "\n".join(f"{s} A G" for s in keep) + "\n")

    # route 1: harmonize all, restrict at export (the pipeline)
    raw_df.to_csv(tmp_path / "all.tsv", sep="\t", index=False)
    run(["python3", BIN / "harmonize_sumstats.py", "--input", tmp_path / "all.tsv", "--id", "t",
         "--out", tmp_path / "all"])
    run(["python3", BIN / "export_sumstats.py", "--input", tmp_path / "all.harmonized.tsv.gz",
         "--format", "lava", "--out", tmp_path / "pipeline.lava.tsv",
         "--snp-list", tmp_path / "list.txt"])

    # route 2: filter the raw file first (manual), then harmonize + export
    raw_df[raw_df.SNP.isin(keep)].to_csv(tmp_path / "sub.tsv", sep="\t", index=False)
    run(["python3", BIN / "harmonize_sumstats.py", "--input", tmp_path / "sub.tsv", "--id", "t",
         "--out", tmp_path / "sub"])
    run(["python3", BIN / "export_sumstats.py", "--input", tmp_path / "sub.harmonized.tsv.gz",
         "--format", "lava", "--out", tmp_path / "manual.lava.tsv"])

    a = (tmp_path / "pipeline.lava.tsv").read_text()
    b = (tmp_path / "manual.lava.tsv").read_text()
    assert a == b
    assert a.count("\n") == 1 + 3          # header + rs2, rs3, rs7


def test_snp_list_with_no_overlap_fails_clearly(tmp_path):
    h = pd.DataFrame({c: ["NA"] for c in hs.HARMONIZED_COLUMNS})
    h.loc[0, ["SNP", "A1", "A2", "Z", "N", "P"]] = ["rs1", "A", "G", "1", "10", "0.3"]
    h.to_csv(tmp_path / "h.tsv.gz", sep="\t", index=False, compression="gzip")
    (tmp_path / "list.txt").write_text("rs999\n")
    with pytest.raises(AssertionError, match="no variants left after --snp-list"):
        run(["python3", BIN / "export_sumstats.py", "--input", tmp_path / "h.tsv.gz",
             "--format", "lava", "--out", tmp_path / "x.tsv", "--snp-list", tmp_path / "list.txt"])
