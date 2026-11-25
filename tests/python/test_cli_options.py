"""True/false options given on the command line, run the way users run the pipeline.

With Nextflow's default (v2) syntax parser, `--sumher_genomic_control false` on
the command line reaches the pipeline as the text "false", which Groovy treats
as true. nf-test cannot catch this: it needs NXF_SYNTAX_PARSER=v1, under which
the same option arrives as a real false. So this test runs
`nextflow run . -profile test,docker -stub` itself, with the default parser.

It needs nextflow and docker, takes a few minutes (it stages the test data),
and only runs when GENCORR_NEXTFLOW_TESTS=1:

    GENCORR_NEXTFLOW_TESTS=1 pytest tests/python/test_cli_options.py

The save_* checks run LDSC + LAVA only. The SumHer check (LDAK's
--genomic-control) is separate and skipped when the LDAK image is not
available locally (e.g. on ARM64, where LDAK cannot be built).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess

import pytest

from conftest import REPO

pytestmark = pytest.mark.skipif(
    os.environ.get("GENCORR_NEXTFLOW_TESTS") != "1"
    or shutil.which("nextflow") is None or shutil.which("docker") is None,
    reason="set GENCORR_NEXTFLOW_TESTS=1 (needs nextflow and docker)")


def run_stub(tmp_path, name: str, *options: str):
    outdir = tmp_path / name
    env = {k: v for k, v in os.environ.items() if k != "NXF_SYNTAX_PARSER"}   # default parser
    r = subprocess.run(["nextflow", "-q", "run", str(REPO), "-profile", "test,docker", "-stub",
                        "--outdir", str(outdir), "-w", str(tmp_path / f"work_{name}"), *options],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return outdir


def genomic_control(outdir) -> set[str]:
    """--genomic-control values passed to LDAK (the stub log records its options)."""
    logs = list((outdir / "sumher" / "logs").glob("*.log"))
    assert logs
    return {m for f in logs for m in re.findall(r"--genomic-control (\w+)", f.read_text())}


LDAK_IMAGE = "ghcr.io/mauricioroza/gencorr-ldak:6.1"   # as in modules/local/sumher/*/main.nf


def have_image(image: str) -> bool:
    if shutil.which("docker") is None:
        return False
    return subprocess.run(["docker", "image", "inspect", image],
                          capture_output=True).returncode == 0


requires_ldak = pytest.mark.skipif(not have_image(LDAK_IMAGE),
                                   reason=f"{LDAK_IMAGE} not available locally (x86_64 only)")


def test_save_options_false_on_the_command_line(tmp_path):
    out = run_stub(tmp_path, "false", "--methods", "ldsc,lava",
                   "--save_intermediates", "false", "--save_harmonized", "false")
    assert not (out / "lava" / "chunks").exists()
    assert all(f.name.endswith(".qc.json") for f in (out / "harmonized").iterdir())


def test_save_options_true_on_the_command_line(tmp_path):
    out = run_stub(tmp_path, "true", "--methods", "ldsc,lava", "--save_intermediates", "true")
    assert (out / "lava" / "chunks").exists()


def test_numeric_qc_options_on_the_command_line(tmp_path):
    # with the default parser "--min_info 0.9" arrives as the text "0.9"; comparing it
    # with 0 in conf/modules.config stopped every harmonization task (2026-11-25)
    out = run_stub(tmp_path, "qc", "--methods", "ldsc", "--min_info", "0.9", "--min_maf", "0.01")
    qc = json.loads((out / "harmonized" / "bmi.qc.json").read_text())
    assert "--min-info 0.9" in qc["stub_args"] and "--min-maf 0.01" in qc["stub_args"]


@requires_ldak
@pytest.mark.parametrize("value,expected", [("false", "NO"), ("true", "YES")])
def test_sumher_genomic_control_on_the_command_line(tmp_path, value, expected):
    out = run_stub(tmp_path, value, "--methods", "sumher", "--sumher_genomic_control", value)
    assert genomic_control(out) == {expected}
