"""rsID mapping in bin/harmonize_sumstats.py (--rsid-map), for files whose SNP IDs
are not rsIDs (e.g. GRCh38 'chr:pos:ref:alt' IDs) or that have no SNP column."""

import gzip
import json

import pandas as pd
import pytest

from conftest import BIN, run

import harmonize_sumstats as hs


def raw(**cols):
    return pd.DataFrame({k: [None if v is None else str(v) for v in vals] for k, vals in cols.items()})


def harmonize(df, **kw):
    return hs.harmonize(df, hs.map_columns(list(df.columns), {}), **kw)


MAP = pd.DataFrame({"CHR": ["1", "1", "2", "3"], "BP": [100, 200, 300, 400],
                    "SNP": ["rs1", "rs2", "rs3", "rs4"],
                    "A1": ["A", "C", "A", "G"], "A2": ["G", "T", "C", "T"]})


def gwas(snp, chrom, bp, a1, a2):
    n = len(snp)
    return raw(SNP=snp, CHR=chrom, BP=bp, A1=a1, A2=a2, Z=[1.0] * n, N=[1000] * n)


def test_position_and_alleles_give_rsid():
    df = gwas(["1:100:A:G", "1:200:C:T"], ["1", "chr1"], [100, 200], ["A", "C"], ["G", "T"])
    out, qc = harmonize(df, rsid_map=MAP, rsid_map_name="map.tsv", build="GRCh38")
    assert list(out.SNP) == ["rs1", "rs2"]
    assert qc["rsid_mapping"] == {"map": "map.tsv", "build": "GRCh38", "rows_without_rsid": 2,
                                  "mapped": 2, "matched_by": "chr, position and alleles (either order or strand)",
                                  "unmatched": 0}


def test_allele_order_and_strand_do_not_matter():
    # swapped alleles (G/A) and the complementary strand (T/C for A/G) both match rs1
    df = gwas(["x1", "x2"], ["1", "1"], [100, 100], ["G", "T"], ["A", "C"])
    out, _ = harmonize(df.iloc[[0]].reset_index(drop=True), rsid_map=MAP)
    assert list(out.SNP) == ["rs1"]
    out, _ = harmonize(df.iloc[[1]].reset_index(drop=True), rsid_map=MAP)
    assert list(out.SNP) == ["rs1"]


def test_different_alleles_at_same_position_do_not_match():
    df = gwas(["1:100:A:T"], ["1"], [100], ["A"], ["T"])   # rs1 at 1:100 is A/G
    out, qc = harmonize(df, rsid_map=MAP)
    assert len(out) == 0 and qc["removed"]["no_rsid_after_mapping"] == 1


def test_unmatched_rows_are_dropped_and_rsids_kept():
    df = gwas(["rs9", "1:100:A:G", "5:999:A:G"], ["1", "1", "5"], [50, 100, 999],
              ["A", "A", "A"], ["G", "G", "G"])
    out, qc = harmonize(df, rsid_map=MAP)
    assert list(out.SNP) == ["rs9", "rs1"]            # existing rsID untouched
    assert qc["rsid_mapping"]["mapped"] == 1 and qc["rsid_mapping"]["unmatched"] == 1
    assert qc["removed"]["no_rsid_after_mapping"] == 1


def test_no_snp_column_works_with_map():
    df = raw(CHR=["2", "3"], BP=[300, 400], A1=["C", "T"], A2=["A", "G"], Z=[1.0, 2.0], N=[10, 10])
    out, qc = harmonize(df, rsid_map=MAP)
    assert list(out.SNP) == ["rs3", "rs4"]


def test_no_snp_column_without_map_fails():
    df = raw(CHR=["2"], BP=[300], A1=["C"], A2=["A"], Z=[1.0], N=[10])
    with pytest.raises(SystemExit):
        harmonize(df)


def test_map_without_alleles_skips_ambiguous_positions():
    m = pd.DataFrame({"CHR": ["1", "1", "1"], "BP": [100, 200, 200], "SNP": ["rs1", "rs2", "rs22"]})
    df = gwas(["a", "b"], ["1", "1"], [100, 200], ["A", "C"], ["G", "T"])
    out, qc = harmonize(df, rsid_map=m)
    assert list(out.SNP) == ["rs1"]                   # position 200 has two rsIDs
    assert qc["rsid_mapping"]["matched_by"].startswith("chr and position")


def test_warning_when_most_ids_are_not_rsids_and_no_map():
    df = gwas(["1:100:A:G", "1:200:C:T", "rs5"], ["1", "1", "1"], [100, 200, 300],
              ["A", "C", "A"], ["G", "T", "G"])
    _, qc = harmonize(df, build="GRCh38")
    assert "warning_non_rsid" in qc and "GRCh38" in qc["warning_non_rsid"]


def test_no_warning_and_no_new_qc_keys_for_rsid_files():
    # keeps QC reports (and test snapshots) of ordinary rsID files unchanged
    df = gwas(["rs1", "rs2"], ["1", "1"], [100, 200], ["A", "C"], ["G", "T"])
    _, qc = harmonize(df)
    assert "warning_non_rsid" not in qc and "rsid_mapping" not in qc
    assert "no_rsid_after_mapping" not in qc["removed"]


def test_read_rsid_map_bim_and_vcf_like_table(tmp_path):
    bim = tmp_path / "ref.bim"
    bim.write_text("1\trs1\t0\t100\tA\tG\n1\t1:150_A_C\t0\t150\tA\tC\nchr2\trs3\t0\t300\tC\tA\n")
    m = hs.read_rsid_map(bim)
    assert list(m.SNP) == ["rs1", "rs3"] and list(m.CHR) == ["1", "2"]   # non-rs rows dropped

    vcf_like = tmp_path / "map.tsv.gz"
    with gzip.open(vcf_like, "wt") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\n1\t100\trs1\tG\tA\n")
    m = hs.read_rsid_map(vcf_like)
    assert m.iloc[0].to_dict() == {"CHR": "1", "BP": 100, "SNP": "rs1", "A1": "A", "A2": "G"}


def test_cli_with_bim_map(tmp_path):
    (tmp_path / "ref.bim").write_text("1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tC\tT\n")
    (tmp_path / "g.tsv").write_text(
        "variant_id\tchromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tp_value\tn\n"
        "1_100_G_A\t1\t100\tA\tG\t0.1\t0.01\t1000\n"
        "1_200_T_C\t1\t200\tC\tT\t-0.2\t0.001\t1000\n")
    run(["python3", BIN / "harmonize_sumstats.py", "--input", tmp_path / "g.tsv", "--id", "t",
         "--out", tmp_path / "t", "--rsid-map", tmp_path / "ref.bim", "--build", "GRCh38"])
    out = pd.read_csv(tmp_path / "t.harmonized.tsv.gz", sep="\t")
    assert list(out.SNP) == ["rs1", "rs2"]
    qc = json.loads((tmp_path / "t.qc.json").read_text())
    assert qc["rsid_mapping"]["mapped"] == 2 and qc["rsid_mapping"]["build"] == "GRCh38"
