"""bin/fetch_gwascatalog.py and GWAS Catalog input in harmonize_sumstats.py.

The metadata (*-meta.yaml) and folder listings in fixtures/gwascatalog are real
GWAS Catalog files (downloaded 2026-11-10). Downloads are tested against a local
folder laid out like the FTP, so no network is needed.
"""

import gzip
import hashlib
import json
import re

import pandas as pd
import pytest

from conftest import BIN, FIXTURES, run

import fetch_gwascatalog as fg

GC = FIXTURES / "gwascatalog"


# ---------------------------------------------------------------- locating --
@pytest.mark.parametrize("acc, expected", [
    ("GCST90000016", "GCST90000001-GCST90001000/GCST90000016"),
    ("GCST000028", "GCST000001-GCST001000/GCST000028"),
    ("GCST006465", "GCST006001-GCST007000/GCST006465"),
    ("GCST001000", "GCST000001-GCST001000/GCST001000"),
    ("GCST90475810", "GCST90475001-GCST90476000/GCST90475810"),
])
def test_study_dir_matches_ftp_layout(acc, expected):
    assert fg.study_dir(acc) == expected


def test_not_an_accession_is_rejected():
    with pytest.raises(SystemExit, match="not a GWAS Catalog accession"):
        fg.study_dir("ieu-a-25")


def test_newer_file_is_preferred_when_a_study_has_two():
    html = (GC / "GCST006465_harmonised_index.html").read_text()
    names = sorted(set(re.findall(r'href="([^"/?]+)"', html)))
    assert sum(n.endswith(".h.tsv.gz") for n in names) == 2       # real case on the FTP
    assert fg.pick_file("GCST006465", names) == "GCST006465.h.tsv.gz"


def test_single_file_is_used_whatever_its_name():
    assert fg.pick_file("GCST000028", ["17463246-GCST000028-EFO_0001360.h.tsv.gz", "readme.txt"]) \
        == "17463246-GCST000028-EFO_0001360.h.tsv.gz"


# ---------------------------------------------------------------- metadata --
def test_metadata_unindented_sample_list_multi_ancestry():
    m = fg.parse_meta_yaml((GC / "GCST90000016.h.tsv.gz-meta.yaml").read_text())
    assert m["data_file_name"] == "GCST90000016.h.tsv.gz"
    assert m["data_file_md5sum"] == "71375cc8bf0acfc57ea8854cd9061e5e"
    assert m["genome_assembly"] == "GRCh38"
    assert m["trait_description"] == ["Migraine"]
    assert len(m["samples"]) == 7
    assert m["samples"][0] == {"ancestry": ["European"], "size": 513266, "case_control": True}
    n, reason, warnings = fg.total_n(m["samples"], "EUR")
    assert n is None and "multi-ancestry" in reason and "European 513266" in reason
    assert warnings == []


def test_metadata_indented_sample_list_single_ancestry():
    m = fg.parse_meta_yaml((GC / "17463246-GCST000028-EFO_0001360.h.tsv.gz-meta.yaml").read_text())
    assert m["data_file_md5sum"] == "0a873f1c0427f2ef31df29851e7ba8ae"
    assert m["samples"] == [{"ancestry": ["European"], "size": 2931, "case_control": True}]
    assert fg.total_n(m["samples"], "EUR") == (2931, None, [])


def test_metadata_value_continued_on_next_line():
    m = fg.parse_meta_yaml((GC / "GCST90475810.h.tsv.gz-meta.yaml").read_text())
    assert m["file_type"] == "GWAS-SSF v1.0"
    assert m["samples"] == [{"ancestry": ["European"], "size": 447460, "case_control": True}]


def test_several_sample_groups_give_no_n():
    # real case (GCST90002409): the same European sample listed twice
    samples = [{"ancestry": ["European"], "size": 39620, "case_control": False}] * 2
    n, reason, _ = fg.total_n(samples, "EUR")
    assert n is None and "2 sample groups" in reason


def test_not_reported_ancestry_gives_no_n():
    # real case (GCST004131): categories NR and European in one entry
    n, reason, warnings = fg.total_n([{"ancestry": ["NR", "European"], "size": 59957,
                                      "case_control": True}], "EUR")
    assert n is None and "NR/European 59957" in reason and warnings == []


def test_ancestry_mismatch_warns_and_gives_no_n():
    samples = [{"ancestry": ["European"], "size": 1000, "case_control": False}]
    n, reason, warnings = fg.total_n(samples, "EAS")
    assert n is None and reason
    assert "does not include East Asian" in warnings[0]


def test_african_categories_count_as_afr():
    # African American and Sub-Saharan African groups are both AFR in the samplesheet
    samples = [{"ancestry": ["African American or Afro-Caribbean"], "size": 900, "case_control": True}]
    assert fg.total_n(samples, "AFR") == (900, None, [])
    n, _, warnings = fg.total_n(samples, "SAS")
    assert n is None and "does not include South Asian" in warnings[0]


def test_multi_ancestry_never_takes_catalog_n():
    samples = [{"ancestry": ["European"], "size": 1000, "case_control": False}]
    n, reason, warnings = fg.total_n(samples, "MULTI")
    assert n is None and "multi-ancestry" in reason and warnings == []


# ------------------------------------------------- download (local mirror) --
def make_study(root, acc, rows, yaml_md5=None, samples="  - sample_ancestry_category:\n"
               "      - European\n    sample_size: 5000\n    case_control_study: false\n"):
    d = root / fg.study_dir(acc) / "harmonised"
    d.mkdir(parents=True)
    data = d / f"{acc}.h.tsv.gz"
    with gzip.open(data, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    md5 = hashlib.md5(data.read_bytes()).hexdigest()
    (d / f"{acc}.h.tsv.gz-meta.yaml").write_text(
        f"gwas_id: {acc}\ngenome_assembly: GRCh38\ntrait_description:\n  - Test trait\n"
        f"samples:\n{samples}data_file_name: {acc}.h.tsv.gz\n"
        f"data_file_md5sum: {yaml_md5 or md5}\n")
    return md5


NEW_FORMAT = ["chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error"
              "\teffect_allele_frequency\tp_value\tvariant_id\trsid",
              "1\t1000\tA\tG\t0.1\t0.05\t0.3\t0.045\t1_1000_G_A\trs1",
              "1\t2000\tC\tT\t-0.2\t0.05\t0.4\t0.0001\t1_2000_T_C\trs2",
              "1\t3000\tA\tC\t0.05\t0.05\t0.2\t0.32\t1_3000_C_A\trs3"]


def test_download_checks_md5_and_writes_source_info(tmp_path):
    mirror = tmp_path / "ftp"
    md5 = make_study(mirror, "GCST90000123", NEW_FORMAT)
    work = tmp_path / "work"
    work.mkdir()
    run(["python3", BIN / "fetch_gwascatalog.py", "--accession", "GCST90000123", "--id", "trait",
         "--ancestry", "EUR", "--base", mirror], cwd=work)
    info = json.loads((work / "trait.gwascatalog.json").read_text())
    assert (work / "GCST90000123.h.tsv.gz").is_file()
    assert info["md5"] == md5 and info["md5_checked"]
    assert info["n"] == 5000 and info["n_reason"] is None
    assert info["trait"] == "Test trait" and info["genome_assembly"] == "GRCh38"


def test_md5_mismatch_fails(tmp_path):
    mirror = tmp_path / "ftp"
    make_study(mirror, "GCST90000123", NEW_FORMAT, yaml_md5="0" * 32)
    with pytest.raises(AssertionError, match="differs from the GWAS Catalog"):
        run(["python3", BIN / "fetch_gwascatalog.py", "--accession", "GCST90000123", "--id", "t",
             "--base", mirror], cwd=tmp_path)


def test_unknown_accession_fails_clearly(tmp_path):
    (tmp_path / "ftp").mkdir()
    with pytest.raises(AssertionError, match="no harmonised summary statistics"):
        run(["python3", BIN / "fetch_gwascatalog.py", "--accession", "GCST99999999", "--id", "t",
             "--base", tmp_path / "ftp"], cwd=tmp_path)


def test_harmonised_list_fallback_matches_whole_accession(tmp_path):
    # study not in its expected folder; harmonised_list.txt also lists a study whose
    # accession merely contains the one asked for (a substring match would take that one)
    mirror = tmp_path / "ftp"
    other = mirror / "elsewhere/GCST0001234/harmonised"
    wanted = mirror / "elsewhere/GCST000123/harmonised"
    for d, val in ((other, "0.9"), (wanted, "0.1")):
        d.mkdir(parents=True)
        with gzip.open(d / f"{d.parent.name}.h.tsv.gz", "wt") as fh:
            fh.write(NEW_FORMAT[0] + "\n" + NEW_FORMAT[1].replace("\t0.1\t", f"\t{val}\t") + "\n")
    (mirror / "harmonised_list.txt").write_text(
        "./elsewhere/GCST0001234/harmonised/GCST0001234.h.tsv.gz\n"
        "./elsewhere/GCST000123/harmonised/GCST000123.h.tsv.gz\n")
    run(["python3", BIN / "fetch_gwascatalog.py", "--accession", "GCST000123", "--id", "t",
         "--base", mirror], cwd=tmp_path)
    info = json.loads((tmp_path / "t.gwascatalog.json").read_text())
    assert info["url"].endswith("elsewhere/GCST000123/harmonised/GCST000123.h.tsv.gz")
    assert any("no metadata file" in w for w in info["warnings"])


# --------------------------------------------- harmonizing a Catalog file --
def harmonize(tmp_path, rows, info, *extra):
    src = tmp_path / "in.h.tsv.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    (tmp_path / "info.json").write_text(json.dumps(info))
    run(["python3", BIN / "harmonize_sumstats.py", "--input", src, "--id", "t", "--out", tmp_path / "t",
         "--source-info", tmp_path / "info.json", *extra])
    qc = json.loads((tmp_path / "t.qc.json").read_text())
    return pd.read_csv(tmp_path / "t.harmonized.tsv.gz", sep="\t"), qc


INFO = {"source": "gwas_catalog", "accession": "GCST90000123", "n": 5000, "n_reason": None,
        "warnings": [], "samples": [], "md5": "x", "md5_checked": True}


def test_n_from_catalog_metadata_and_provenance_in_qc(tmp_path):
    out, qc = harmonize(tmp_path, NEW_FORMAT, INFO)
    assert (out["N"] == 5000).all()
    assert qc["n_source"] == "gwas_catalog_metadata"
    assert qc["source"]["accession"] == "GCST90000123"
    assert list(out["SNP"]) == ["rs1", "rs2", "rs3"]                  # harmonised rsid, not variant_id


def test_samplesheet_n_wins_over_catalog_metadata(tmp_path):
    out, qc = harmonize(tmp_path, NEW_FORMAT, INFO, "--n", "1234")
    assert (out["N"] == 1234).all() and qc["n_source"] == "samplesheet"


def test_no_usable_n_explains_why(tmp_path):
    info = {**INFO, "n": None, "n_reason": "multi-ancestry or other-ancestry study (...)"}
    with pytest.raises(AssertionError, match="GWAS Catalog N cannot be used: multi-ancestry"):
        harmonize(tmp_path, NEW_FORMAT, info)


OLD_FORMAT = [  # submitted columns next to harmonised hm_ columns; here the
                # harmonisation flipped the alleles of rs2 (effect sign flips too)
    "hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta"
    "\thm_effect_allele_frequency\thm_code\tvariant_id\tchromosome\tbase_pair_location"
    "\tother_allele\teffect_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value",
    "1_1000_G_A\trs1\t1\t1000\tG\tA\t0.1\t0.3\t10\trs1\t1\t900\tG\tA\t0.1\t0.05\t0.3\t0.045",
    "1_2000_C_T\trs2\t1\t2000\tC\tT\t0.2\t0.6\t11\trs2\t1\t1900\tT\tC\t-0.2\t0.05\t0.4\t0.0001",
]


def test_old_format_uses_the_harmonised_columns_together(tmp_path):
    out, qc = harmonize(tmp_path, OLD_FORMAT, INFO)
    assert qc["column_mapping"]["A1"] == "hm_effect_allele"
    assert qc["column_mapping"]["BETA"] == "hm_beta"
    r2 = out.set_index("SNP").loc["rs2"]
    assert (r2["A1"], r2["A2"]) == ("T", "C") and r2["Z"] > 0          # same association either way


def test_colmap_overrides_the_harmonised_columns(tmp_path):
    _, qc = harmonize(tmp_path, OLD_FORMAT, INFO, "--colmap",
                      "A1=effect_allele,A2=other_allele,BETA=beta,EAF=effect_allele_frequency")
    assert qc["column_mapping"]["A1"] == "effect_allele"


# ------------------------------------------------------------ HTTP download --
class _RangeServer:
    """Local HTTP server for one file: serves byte ranges (or not) and can drop
    every connection after `drop_after` bytes, like a slow, flaky FTP mirror."""

    def __init__(self, data: bytes, ranges: bool = True, drop_after: int | None = None):
        import http.server
        import threading

        outer = self
        self.data, self.ranges, self.drop_after = data, ranges, drop_after
        self.range_requests = 0

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _span(self):
                rng = self.headers.get("Range") if outer.ranges else None
                if not rng:
                    return 0, len(outer.data) - 1, False
                outer.range_requests += 1
                a, b = rng.split("=")[1].split("-")
                return int(a), (int(b) if b else len(outer.data) - 1), True

            def do_HEAD(self):
                self.send_response(200)
                self.send_header("Content-Length", str(len(outer.data)))
                if outer.ranges:
                    self.send_header("Accept-Ranges", "bytes")
                self.end_headers()

            def do_GET(self):
                a, b, partial = self._span()
                body = outer.data[a:b + 1]
                self.send_response(206 if partial else 200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if outer.drop_after is not None:
                    body = body[:outer.drop_after]  # connection closes early
                self.wfile.write(body)
                self.close_connection = True

        self.httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def close(self):
        self.httpd.shutdown()


@pytest.fixture
def payload():
    import random
    return random.Random(1).randbytes(300_000)


@pytest.mark.parametrize("connections", [1, 4])
def test_http_download_in_ranges_matches_the_file(tmp_path, monkeypatch, payload, connections):
    monkeypatch.setattr(fg, "MIN_SPLIT", 1000)
    srv = _RangeServer(payload)
    try:
        md5 = fg.download(srv.url, "f.gz", tmp_path / "f.gz", connections)
    finally:
        srv.close()
    assert (tmp_path / "f.gz").read_bytes() == payload
    assert md5 == hashlib.md5(payload).hexdigest()
    if connections == 4:
        assert srv.range_requests == 4


def test_http_download_continues_after_dropped_connections(tmp_path, monkeypatch, payload):
    monkeypatch.setattr(fg, "MIN_SPLIT", 1000)
    monkeypatch.setattr(fg.time, "sleep", lambda _s: None)
    srv = _RangeServer(payload, drop_after=20_000)  # every connection stops after 20 kB
    try:
        md5 = fg.download(srv.url, "f.gz", tmp_path / "f.gz", 2)
    finally:
        srv.close()
    assert md5 == hashlib.md5(payload).hexdigest()
    assert srv.range_requests >= len(payload) // 20_000  # resumed, never restarted


def test_http_download_without_range_support(tmp_path, monkeypatch, payload):
    monkeypatch.setattr(fg, "MIN_SPLIT", 1000)
    srv = _RangeServer(payload, ranges=False)
    try:
        md5 = fg.download(srv.url, "f.gz", tmp_path / "f.gz", 4)
    finally:
        srv.close()
    assert md5 == hashlib.md5(payload).hexdigest() and srv.range_requests == 0


def test_http_download_gives_up_when_no_data_arrives(tmp_path, monkeypatch, payload):
    monkeypatch.setattr(fg.time, "sleep", lambda _s: None)
    srv = _RangeServer(payload, drop_after=0)
    try:
        with pytest.raises(SystemExit):
            fg.download(srv.url, "f.gz", tmp_path / "f.gz", 1)
    finally:
        srv.close()
