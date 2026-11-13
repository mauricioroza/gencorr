#!/usr/bin/env python3
"""
compare_ldsc_impl.py - compare two gencorr runs that differ only in --ldsc_impl
(LDSC 1.0.1 vs the CBIIT/ldsc Python 3 port), output by output.

Compared (both runs need --save_intermediates true for the munged files):
  ldsc/gencorr.ldsc.tsv           every column, incl. the JSON `extra` fields
  ldsc/gencorr.ldsc_h2.tsv        every column
  ldsc/gencorr.sample_overlap.txt every cell
  ldsc/munged/<trait>.sumstats.gz SNP, A1, A2, Z, N row by row
  ldsc/munged/<trait>.munge.log   SNP counts reported by munge

Numbers are compared as printed by LDSC (4 decimals, p in %.4g), so "identical"
means identical output, not just within tolerance.

    compare_ldsc_impl.py --a run_original --b run_py3 [--report prefix]

Exit status 0 when everything is identical, 1 otherwise.
"""

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def expand_extra(df):
    if "extra" not in df.columns:
        return df
    extra = pd.DataFrame([json.loads(x) if x else {} for x in df["extra"]], index=df.index)
    extra = extra.astype(str).add_prefix("extra.")
    return pd.concat([df.drop(columns="extra"), extra], axis=1)


def compare_frames(a, b, keys, label, out):
    """Return number of differing cells; append lines to out."""
    a = expand_extra(a).set_index(keys).sort_index()
    b = expand_extra(b).set_index(keys).sort_index()
    if not a.index.equals(b.index):
        out.append(f"{label}: different rows: only A {sorted(set(a.index) - set(b.index))}, "
                   f"only B {sorted(set(b.index) - set(a.index))}")
        return 1
    cols = sorted(set(a.columns) | set(b.columns))
    n_diff, max_abs = 0, 0.0
    for c in cols:
        if c not in a.columns or c not in b.columns:
            out.append(f"{label}: column {c} only in {'A' if c in a.columns else 'B'}")
            n_diff += 1
            continue
        for idx in a.index:
            x, y = a.at[idx, c], b.at[idx, c]
            if x == y:
                continue
            n_diff += 1
            try:
                d = abs(float(x) - float(y))
                max_abs = max(max_abs, d)
            except ValueError:
                d = None
            out.append(f"{label}: {idx} {c}: A={x} B={y}" + (f" (|diff| {d:.3g})" if d is not None else ""))
    out.append(f"{label}: {len(a)} rows x {len(cols)} columns, {n_diff} cells differ"
               + (f", max |diff| {max_abs:.3g}" if n_diff else ""))
    return n_diff


def read_munged(path):
    return pd.read_csv(path, sep="\t", dtype={"SNP": str, "A1": str, "A2": str})


def munge_counts(log):
    text = Path(log).read_text()
    pat = re.compile(r"^(Read \d+ SNPs.*|Removed \d+ SNPs.*|Writing summary statistics for \d+ SNPs.*"
                     r"|After merging with --merge-alleles.*|\d+ SNPs remain.*)$", re.M)
    # Python 2 prints floats with 12 significant digits, Python 3 with 17
    # (N < 71533.3333333 vs 71533.33333333333): compare at 10 digits.
    num = re.compile(r"\d+\.\d+")
    return [num.sub(lambda m: f"{float(m.group(0)):.10g}", m.group(1).strip()) for m in pat.finditer(text)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, type=Path, help="gencorr outdir, --ldsc_impl original")
    ap.add_argument("--b", required=True, type=Path, help="gencorr outdir, --ldsc_impl py3")
    ap.add_argument("--report", help="write the report to <prefix>.txt")
    args = ap.parse_args()

    out = [f"A: {args.a}", f"B: {args.b}", ""]
    n = 0
    n += compare_frames(read_tsv(args.a / "ldsc/gencorr.ldsc.tsv"), read_tsv(args.b / "ldsc/gencorr.ldsc.tsv"),
                        ["trait1", "trait2"], "ldsc.tsv", out)
    n += compare_frames(read_tsv(args.a / "ldsc/gencorr.ldsc_h2.tsv"), read_tsv(args.b / "ldsc/gencorr.ldsc_h2.tsv"),
                        [read_tsv(args.a / "ldsc/gencorr.ldsc_h2.tsv").columns[0]], "ldsc_h2.tsv", out)

    oa = pd.read_csv(args.a / "ldsc/gencorr.sample_overlap.txt", sep=r"\s+", index_col=0)
    ob = pd.read_csv(args.b / "ldsc/gencorr.sample_overlap.txt", sep=r"\s+", index_col=0)
    same = oa.shape == ob.shape and np.array_equal(oa.values, ob.values)
    out.append(f"sample_overlap: {'identical' if same else 'DIFFERENT'} ({oa.shape[0]}x{oa.shape[1]})")
    n += 0 if same else 1

    for fa in sorted((args.a / "ldsc/munged").glob("*.sumstats.gz")):
        fb = args.b / "ldsc/munged" / fa.name
        if not fb.exists():
            out.append(f"munged {fa.name}: missing in B")
            n += 1
            continue
        ma, mb = read_munged(fa), read_munged(fb)
        if ma.shape != mb.shape or not (ma[["SNP", "A1", "A2"]].fillna("") == mb[["SNP", "A1", "A2"]].fillna("")).all().all():
            out.append(f"munged {fa.name}: SNP/allele rows differ ({len(ma)} vs {len(mb)} rows)")
            n += 1
            continue
        dz = np.nanmax(np.abs(ma.Z.values - mb.Z.values)) if len(ma) else 0.0
        dn = np.nanmax(np.abs(ma.N.values - mb.N.values)) if len(ma) else 0.0
        na_same = (ma.Z.isna() == mb.Z.isna()).all()
        ok = dz == 0 and dn == 0 and na_same
        n += 0 if ok else 1
        out.append(f"munged {fa.name}: {len(ma)} rows, {ma.Z.notna().sum()} with Z; "
                   f"max |dZ| {dz:.3g}, max |dN| {dn:.3g}, NA pattern {'same' if na_same else 'DIFFERENT'}")
        la, lb = fa.with_name(fa.name.replace(".sumstats.gz", ".munge.log")), fb.with_name(fb.name.replace(".sumstats.gz", ".munge.log"))
        if la.exists() and lb.exists():
            ca, cb = munge_counts(la), munge_counts(lb)
            if ca != cb:
                n += 1
                out.append(f"  munge log counts differ:\n    A: {ca}\n    B: {cb}")

    out.append("")
    out.append("RESULT: " + ("IDENTICAL" if n == 0 else f"{n} difference(s)"))
    text = "\n".join(out)
    print(text)
    if args.report:
        Path(args.report + ".txt").write_text(text + "\n")
    return 0 if n == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
