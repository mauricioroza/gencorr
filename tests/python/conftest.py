"""
Shared fixtures for the Python test suite.

Test layers (run all with `pytest tests/python`):

  unit         pure-Python checks of bin/ scripts on small synthetic data
               (always run, need only pandas/scipy)
  ldsc         oracle tests against the real LDSC 1.0.1 (need ldsc.py and
               munge_sumstats.py on PATH; skipped otherwise)
  lava         oracle tests against the real LAVA R package (need Rscript with
               LAVA installed; skipped otherwise)

The LDSC/LAVA tests download the LAVA vignette data (pinned commit) once into
$GENCORR_TEST_CACHE (default ~/.cache/gencorr-test).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(BIN))

LAVA_COMMIT = "e729a245f7b6923967a96804fbf5246eadf2d6c6"
VIGNETTE = f"https://raw.githubusercontent.com/josefin-werme/LAVA/{LAVA_COMMIT}/vignettes/data"


def have_ldsc() -> bool:
    return shutil.which("ldsc.py") is not None and shutil.which("munge_sumstats.py") is not None


def have_lava() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('LAVA', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


requires_ldsc = pytest.mark.skipif(not have_ldsc(), reason="ldsc.py / munge_sumstats.py not on PATH")
requires_lava = pytest.mark.skipif(not have_lava(), reason="Rscript with the LAVA package not available")


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Run a command; on failure show its output in the pytest report."""
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n"
                             f"--- stdout ---\n{r.stdout[-4000:]}\n--- stderr ---\n{r.stderr[-4000:]}")
    return r


@pytest.fixture(scope="session")
def vignette() -> Path:
    """LAVA vignette data (sumstats, 1000G test reference, loci), cached."""
    cache = Path(os.environ.get("GENCORR_TEST_CACHE", Path.home() / ".cache" / "gencorr-test"))
    cache = cache / LAVA_COMMIT[:12]
    cache.mkdir(parents=True, exist_ok=True)
    names = ["bmi.sumstats.txt", "neuro.sumstats.txt", "depression.sumstats.txt",
             "rheuma.sumstats.txt", "g1000_test.bed", "g1000_test.bim", "g1000_test.fam",
             "test.loci"]
    for name in names:
        target = cache / name
        if not target.exists():
            tmp = target.with_suffix(target.suffix + ".part")
            urllib.request.urlretrieve(f"{VIGNETTE}/{name}", tmp)
            tmp.rename(target)
    return cache


@pytest.fixture(scope="session")
def test_ld() -> Path:
    return REPO / "tests" / "data" / "ldsc"
