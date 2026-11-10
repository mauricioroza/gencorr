"""True/false options given on the command line, run the way users run the pipeline.

With Nextflow's default (v2) syntax parser, `--sumher_genomic_control false` on
the command line reaches the pipeline as the text "false", which Groovy treats
as true. nf-test cannot catch this: it needs NXF_SYNTAX_PARSER=v1, under which
the same option arrives as a real false. So this test runs
`nextflow run . -profile test,docker -stub` itself, with the default parser.

It needs nextflow and docker, takes a few minutes (it stages the test data),
and only runs when GENCORR_NEXTFLOW_TESTS=1:

    GENCORR_NEXTFLOW_TESTS=1 pytest tests/python/test_cli_options.py
"""

from __future__ import annotations

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


def test_false_on_the_command_line_means_false(tmp_path):
    out = run_stub(tmp_path, "false", "--sumher_genomic_control", "false",
                   "--save_intermediates", "false", "--save_harmonized", "false")
    assert genomic_control(out) == {"NO"}
    assert not (out / "lava" / "chunks").exists()
    assert all(f.name.endswith(".qc.json") for f in (out / "harmonized").iterdir())


def test_true_on_the_command_line_means_true(tmp_path):
    out = run_stub(tmp_path, "true", "--sumher_genomic_control", "true", "--save_intermediates", "true")
    assert genomic_control(out) == {"YES"}
    assert (out / "lava" / "chunks").exists()
