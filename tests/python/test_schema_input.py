"""Samplesheet patterns in assets/schema_input.json (checked here with Python's
re; the patterns use no syntax where it differs from nf-schema's)."""

import csv
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PROPS = json.loads((ROOT / "assets" / "schema_input.json").read_text())["items"]["properties"]


def matches(field, value):
    return re.search(PROPS[field]["pattern"], value) is not None


@pytest.mark.parametrize("value,ok", [
    ("bmi", True), ("t2d_ukb", True), ("a_b_c", True), ("x.v2-1", True), ("A1", True),
    ("a__b", False),   # '__' joins trait IDs: a__b + c and a + b__c would both give a__b__c
    ("a_", False),     # a_ + b and a + _b... would give a___b
    ("1abc", False), ("_a", False), ("a b", False), ("a/b", False), ("", False),
])
def test_id_pattern(value, ok):
    assert matches("id", value) is ok


@pytest.mark.parametrize("value,ok", [
    ("Z=STAT", True), ("Z=STAT,N=NMISS", True), ("A1=Allele2,A2=Allele1", True),
    ("Z=it's", False),        # would close the --colmap '...' shell quote
    ("Z=a b", False), ("Z=", False), ("z=STAT", False), ("Z=STAT,", False),
])
def test_colmap_pattern(value, ok):
    assert matches("colmap", value) is ok


def test_committed_samplesheets_are_valid():
    sheets = [p for p in ROOT.glob("tests/**/*.csv") if "samplesheet" in p.name]
    sheets += list(ROOT.glob("assets/samplesheet*.csv"))
    assert sheets
    for sheet in sheets:
        with open(sheet, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("id"):
                    assert matches("id", row["id"]), (sheet, row["id"])
                if row.get("colmap"):
                    assert matches("colmap", row["colmap"]), (sheet, row["colmap"])
