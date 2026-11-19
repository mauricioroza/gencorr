"""
bin/gwama.py - N-weighted multivariate GWAMA.

Unit tests on small synthetic data (always run), plus an oracle test against
the authors' R function (github.com/baselmans/multivariate_GWAMA) on their
example data, both fetched at a pinned commit at test time (the repository has
no licence, so nothing from it is copied here). The oracle test needs Rscript
with data.table (e.g. inside ghcr.io/mauricioroza/gencorr-tools) and is
skipped otherwise.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import BIN, run

sys.path.insert(0, str(BIN))
import gwama  # noqa: E402
from gencorr_common import HARMONIZED_COLUMNS, RESULT_COLUMNS, write_tsv  # noqa: E402

SCRIPT = BIN / "gwama.py"
GWAMA_COMMIT = "d55b6171e35cec19c875cb123bcdbb6559fd69af"
GWAMA_RAW = f"https://raw.githubusercontent.com/baselmans/multivariate_GWAMA/{GWAMA_COMMIT}/Test_Data"


def harmonized(path: Path, rows: list[tuple]) -> Path:
    """rows: (SNP, CHR, BP, A1, A2, Z, N, EAF)"""
    df = pd.DataFrame(rows, columns=["SNP", "CHR", "BP", "A1", "A2", "Z", "N", "EAF"])
    df["BETA"] = np.nan
    df["SE"] = np.nan
    df["P"] = np.nan
    df["INFO"] = np.nan
    df[HARMONIZED_COLUMNS].to_csv(path, sep="\t", index=False, na_rep="NA",
                                  compression={"method": "gzip", "mtime": 0})
    return path


def read_out(prefix: Path) -> pd.DataFrame:
    return pd.read_csv(f"{prefix}.harmonized.tsv.gz", sep="\t", dtype={"CHR": str})


def write_cti(path: Path, traits: list[str], mat) -> Path:
    pd.DataFrame(mat, index=traits, columns=traits).to_csv(path, sep="\t")
    return path


def ldsc_tables(tmp: Path, h2: dict, pairs: list[tuple]) -> tuple[Path, Path]:
    """h2: trait -> (h2_obs, intercept); pairs: (t1, t2, rg, gcov_int)"""
    h2_df = pd.DataFrame([{"trait": t, "h2_obs": v[0], "intercept": v[1], "status": "ok"}
                          for t, v in h2.items()])
    rows = []
    for t1, t2, rg, gcov in pairs:
        row = {c: np.nan for c in RESULT_COLUMNS}
        row.update(method="ldsc", scope="global", trait1=t1, trait2=t2, rg=rg, status="ok",
                   extra=json.dumps({"gcov_int": gcov}))
        rows.append(row)
    write_tsv(h2_df, tmp / "gencorr.ldsc_h2.tsv")
    write_tsv(pd.DataFrame(rows, columns=RESULT_COLUMNS), tmp / "gencorr.ldsc.tsv")
    return tmp / "gencorr.ldsc_h2.tsv", tmp / "gencorr.ldsc.tsv"


@pytest.fixture
def two_traits(tmp_path):
    a = harmonized(tmp_path / "a.tsv.gz", [
        ("rs1", "1", 100, "A", "G", 2.0, 1000, 0.2),
        ("rs2", "1", 200, "C", "T", -1.0, 1000, 0.5),
        ("rs3", "2", 300, "A", "C", 3.0, 1000, 0.3),
    ])
    b = harmonized(tmp_path / "b.tsv.gz", [
        ("rs1", "1", 100, "A", "G", 1.0, 3000, 0.4),
        ("rs2", "1", 200, "T", "C", 1.5, 3000, 0.6),   # alleles swapped
        ("rs4", "3", 400, "G", "T", 4.0, 3000, 0.1),   # only in b
    ])
    return a, b


def cli(tmp_path, files, extra) -> Path:
    out = tmp_path / "combo"
    run([sys.executable, SCRIPT, "--id", "combo", "--traits", "a", "b",
         "--sumstats", *files, "--out", out, *extra])
    return out


# --------------------------------------------------------------------------- #
# Formula
# --------------------------------------------------------------------------- #
def test_independent_traits_equal_weights(tmp_path, two_traits):
    """C = I and equal N*h2: Z = (Z1 + Z2) / sqrt(2), N = N1 + N2."""
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], np.eye(2))
    out = read_out(cli(tmp_path, two_traits, ["--cti", cti, "--h2", "0.3", "0.1"])).set_index("SNP")
    # w = N*h2 = 300 for both traits
    assert out.at["rs1", "Z"] == pytest.approx((2.0 + 1.0) / np.sqrt(2))
    assert out.at["rs1", "N"] == 4000
    # rs2: b's alleles are swapped -> its Z is negated before combining
    assert out.at["rs2", "Z"] == pytest.approx((-1.0 - 1.5) / np.sqrt(2))
    assert (out.at["rs2", "A1"], out.at["rs2", "A2"]) == ("C", "T")


def test_correlated_traits_match_formula(tmp_path, two_traits):
    C = np.array([[1.05, 0.2], [0.2, 0.98]])
    h2 = np.array([0.1, 0.05])
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], C)
    out = read_out(cli(tmp_path, two_traits, ["--cti", cti, "--h2", "0.1", "0.05"])).set_index("SNP")
    N = np.array([1000, 3000])
    s = np.sqrt(N * h2)
    z = np.array([2.0, 1.0])
    assert out.at["rs1", "Z"] == pytest.approx(s @ z / np.sqrt(s @ C @ s))
    assert out.at["rs1", "N"] == pytest.approx(np.sqrt(N) @ np.linalg.inv(C) @ np.sqrt(N))
    eaf = (0.2 * 1000 + 0.4 * 3000) / 4000
    assert out.at["rs1", "EAF"] == pytest.approx(eaf)
    se = 1 / np.sqrt(4000 * 2 * eaf * (1 - eaf))
    assert out.at["rs1", "SE"] == pytest.approx(se)
    assert out.at["rs1", "BETA"] == pytest.approx(out.at["rs1", "Z"] * se)


def test_union_snp_in_one_trait(tmp_path, two_traits):
    """A SNP in one trait only: Z / sqrt(C_ii), N = that trait's N."""
    C = np.array([[1.1, 0.2], [0.2, 1.21]])
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], C)
    out = read_out(cli(tmp_path, two_traits, ["--cti", cti, "--h2", "0.1", "0.1"])).set_index("SNP")
    assert out.at["rs4", "Z"] == pytest.approx(4.0 / 1.1)
    assert out.at["rs4", "N"] == 3000
    assert out.at["rs3", "Z"] == pytest.approx(3.0 / np.sqrt(1.1))
    assert set(out.index) == {"rs1", "rs2", "rs3", "rs4"}


def test_sign_reverses_trait_and_its_covariance(tmp_path, two_traits):
    """Flipping b's sign == giving b negated Z-scores and a negated gcov_int."""
    C = np.array([[1.0, -0.3], [-0.3, 1.0]])
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], C)
    flipped = read_out(cli(tmp_path, two_traits,
                           ["--cti", cti, "--h2", "0.1", "0.1", "--signs", "1", "-1"]))
    a, b = two_traits
    neg = pd.read_csv(b, sep="\t")
    neg["Z"] = -neg["Z"]
    neg.to_csv(tmp_path / "bneg.tsv.gz", sep="\t", index=False, na_rep="NA")
    cti2 = write_cti(tmp_path / "cti2.txt", ["a", "b"], [[1.0, 0.3], [0.3, 1.0]])
    manual = read_out(cli(tmp_path, [a, tmp_path / "bneg.tsv.gz"],
                          ["--cti", cti2, "--h2", "0.1", "0.1"]))
    pd.testing.assert_frame_equal(flipped, manual)


def test_allele_mismatch_dropped_for_that_trait(tmp_path):
    a = harmonized(tmp_path / "a.tsv.gz", [("rs1", "1", 1, "A", "G", 2.0, 100, 0.5)])
    b = harmonized(tmp_path / "b.tsv.gz", [("rs1", "1", 1, "A", "C", 9.0, 100, 0.5)])
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], np.eye(2))
    out = cli(tmp_path, [a, b], ["--cti", cti, "--h2", "0.1", "0.1"])
    df = read_out(out)
    assert df["Z"].tolist() == [pytest.approx(2.0)]
    qc = json.loads(Path(f"{out}.gwama.json").read_text())
    assert qc["traits"][1]["allele_mismatch"] == 1


def test_effective_n_independent_is_sum():
    N = np.array([[100.0, 200.0], [100.0, np.nan]])
    assert gwama.effective_n(N, np.eye(2)).tolist() == [300.0, 100.0]


# --------------------------------------------------------------------------- #
# CTI from LDSC
# --------------------------------------------------------------------------- #
def test_cti_from_ldsc_tables(tmp_path):
    h2p, resp = ldsc_tables(tmp_path, {"a": (0.2, 1.02), "b": (0.1, 0.99), "c": (0.3, 1.0)},
                            [("a", "b", 0.8, 0.10), ("b", "a", 0.8, 0.12), ("a", "c", 0.7, 0.0),
                             ("b", "c", 0.9, 0.05)])
    h2, cti, rg = gwama.cti_from_ldsc(h2p, resp, ["a", "b", "c"])
    assert h2.tolist() == [0.2, 0.1, 0.3]
    assert np.diag(cti).tolist() == [1.02, 0.99, 1.0]
    assert cti[0, 1] == cti[1, 0] == pytest.approx(0.11)    # both orientations averaged
    assert rg[1, 2] == 0.9


def test_ldsc_route_and_low_rg_warning(tmp_path, two_traits):
    h2p, resp = ldsc_tables(tmp_path, {"a": (0.2, 1.0), "b": (0.1, 1.0)}, [("a", "b", 0.3, 0.1)])
    out = tmp_path / "combo"
    r = run([sys.executable, SCRIPT, "--id", "combo", "--traits", "a", "b",
             "--sumstats", *two_traits, "--ldsc-h2", h2p, "--ldsc-results", resp, "--out", out])
    assert "rg = 0.300 < 0.5" in r.stderr
    qc = json.loads(Path(f"{out}.gwama.json").read_text())
    assert qc["cti"]["source"] == "ldsc" and qc["cti"]["matrix"][0][1] == 0.1
    assert len(qc["warnings"]) == 1


def test_negative_rg_suggests_sign(tmp_path, two_traits):
    h2p, resp = ldsc_tables(tmp_path, {"a": (0.2, 1.0), "b": (0.1, 1.0)}, [("a", "b", -0.8, -0.1)])
    out = tmp_path / "combo"
    base = [sys.executable, SCRIPT, "--id", "combo", "--traits", "a", "b", "--sumstats",
            *two_traits, "--ldsc-h2", h2p, "--ldsc-results", resp, "--out", out]
    assert "negative after signs" in run(base).stderr
    r = run(base + ["--signs", "1", "-1"])
    assert "WARNING" not in r.stderr
    qc = json.loads(Path(f"{out}.gwama.json").read_text())
    assert qc["rg"]["matrix"][0][1] == 0.8 and qc["cti"]["matrix"][0][1] == 0.1


@pytest.mark.parametrize("h2, pairs, message", [
    ({"a": (0.2, 1.0), "b": (0.1, 1.0)}, [], "no LDSC cross-trait intercept"),
    ({"a": (0.2, 1.0), "b": (-0.01, 1.0)}, [("a", "b", 0.8, 0.1)], "not positive"),
    ({"a": (0.2, np.nan), "b": (0.1, 1.0)}, [("a", "b", 0.8, 0.1)], "no usable LDSC h2 intercept"),
    ({"a": (0.2, 1.0), "b": (0.1, 1.0)}, [("a", "b", 0.8, 1.5)], "not positive definite"),
])
def test_ldsc_input_errors(tmp_path, two_traits, h2, pairs, message):
    h2p, resp = ldsc_tables(tmp_path, h2, pairs)
    r = subprocess.run([sys.executable, str(SCRIPT), "--id", "combo", "--traits", "a", "b",
                        "--sumstats", *map(str, two_traits), "--ldsc-h2", str(h2p),
                        "--ldsc-results", str(resp), "--out", str(tmp_path / "x")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and message in r.stderr


def test_output_is_byte_reproducible(tmp_path, two_traits):
    cti = write_cti(tmp_path / "cti.txt", ["a", "b"], [[1.0, 0.1], [0.1, 1.0]])
    digests = []
    for d in ("r1", "r2"):
        (tmp_path / d).mkdir()
        out = cli(tmp_path / d, two_traits, ["--cti", cti, "--h2", "0.1", "0.2"])
        digests.append(hashlib.md5(Path(f"{out}.harmonized.tsv.gz").read_bytes()).hexdigest())
        with gzip.open(f"{out}.harmonized.tsv.gz", "rt") as fh:
            assert fh.readline().rstrip("\n").split("\t") == HARMONIZED_COLUMNS
    assert digests[0] == digests[1]


# --------------------------------------------------------------------------- #
# Oracle: the authors' R function on their example data
# --------------------------------------------------------------------------- #
def have_r_datatable() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('data.table', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


@pytest.fixture(scope="session")
def gwama_example() -> Path:
    cache = Path(os.environ.get("GENCORR_TEST_CACHE", Path.home() / ".cache" / "gencorr-test"))
    cache = cache / f"multivariate_GWAMA-{GWAMA_COMMIT[:12]}"
    cache.mkdir(parents=True, exist_ok=True)
    names = ["N_weighted_GWAMA.function.1_2_6.R", "cross_trait_intercept.txt",
             "LS_100K_no23andMe.txt", "PA_100K_no23andMe.txt",
             "NEU_100K_no23andMe.txt", "DEP_100K_no23andMe.txt"]
    for name in names:
        target = cache / name
        if not target.exists():
            tmp = target.with_suffix(target.suffix + ".part")
            urllib.request.urlretrieve(f"{GWAMA_RAW}/{name}", tmp)
            tmp.rename(target)
    return cache


@pytest.mark.skipif(not have_r_datatable(), reason="Rscript with data.table not available")
def test_oracle_against_authors_r_function(tmp_path, gwama_example):
    """Same inputs as the authors' example_script.R (h2 0.0498, 0.0441, 0.0708,
    0.0294 and their CTI file); every SNP's Z, P, N_eff, BETA, SE and EAF must
    agree with the R function's output."""
    d = gwama_example
    traits = ["LS", "PA", "NEU", "DEP"]
    h2 = ["0.0498", "0.0441", "0.0708", "0.0294"]
    r_script = tmp_path / "run.R"
    r_script.write_text(f"""
suppressMessages(library(data.table))
source("{d / 'N_weighted_GWAMA.function.1_2_6.R'}")
CTI <- as.matrix(read.table("{d / 'cross_trait_intercept.txt'}", header=TRUE, row.names=1))
dat <- lapply(c({', '.join(f'"{d / (t + "_100K_no23andMe.txt")}"' for t in traits)}),
              function(f) fread(f, showProgress=FALSE, data.table=FALSE))
multivariate_GWAMA(x=dat, cov_Z=CTI, h2=c({', '.join(h2)}), out="{tmp_path}/",
                   name="oracle", output_gz=FALSE, check_columns=FALSE)
""")
    run(["Rscript", r_script], cwd=tmp_path)
    oracle = pd.read_csv(tmp_path / "oracle.N_weighted_GWAMA.results.txt", sep="\t")

    # the authors' files -> harmonized format (A1 = EA). The R function drops
    # every row with a missing value, including a missing EAF (45 rows in three
    # of these files); gwama.py keeps such SNPs, as the harmonized format allows
    # a missing EAF. Drop them here too so both get the same input.
    files = []
    for t in traits:
        raw = pd.read_csv(d / f"{t}_100K_no23andMe.txt", sep=r"\s+").dropna()
        raw.columns = ["SNP", "CHR", "BP", "A1", "A2", "EAF", "N", "Z", "P"]
        for c in ("BETA", "SE", "INFO"):
            raw[c] = np.nan
        f = tmp_path / f"{t}.harmonized.tsv.gz"
        raw[HARMONIZED_COLUMNS].to_csv(f, sep="\t", index=False, na_rep="NA")
        files.append(f)
    out = tmp_path / "py"
    run([sys.executable, SCRIPT, "--id", "WB", "--traits", *traits, "--sumstats", *files,
         "--cti", d / "cross_trait_intercept.txt", "--h2", *h2, "--out", out])
    py = read_out(out)

    merged = oracle.merge(py, left_on="SNPID", right_on="SNP", suffixes=("_r", "_py"))
    assert len(merged) == len(oracle) == len(py)
    assert (merged["EA"].str.upper() == merged["A1"]).all()
    for r_col, py_col in [("Z_r", "Z_py"), ("N_eff", "N"), ("BETA_r", "BETA_py"),
                          ("SE_r", "SE_py"), ("EAF_r", "EAF_py")]:
        np.testing.assert_allclose(merged[py_col], merged[r_col], rtol=1e-6, err_msg=r_col)
    np.testing.assert_allclose(np.log10(merged["P"]), np.log10(merged["PVAL"]), rtol=1e-6)
