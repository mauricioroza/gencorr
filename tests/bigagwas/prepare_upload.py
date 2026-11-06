#!/usr/bin/env python3
"""
prepare_upload.py - turn the 4 LAVA vignette summary statistics into files
that BigaGWAS accepts for upload, so the SAME files can be run through
BigaGWAS and gencorr and the results compared.

BigaGWAS upload rules (bigagwas/process_input.py, general_check_for_harmo):
  * required: a1, a2, se, pvalue, n  +  beta or odds_ratio
  * plus snp (rsID) or chr + position
  * header names are matched case-insensitively against a fixed list;
    N_analyzed / StdErrLogOR / SNPID_UKB-style names are not all recognised.

What each file needs:
  bmi        SNPID_UKB -> SNP, NMISS -> N; drop MAF (BigaGWAS would read it as
             EAF), STAT, INFO_UKB
  neuro      N_analyzed -> N; drop Z, EAF, INFO so both tools derive Z from P
  depression StdErrLogOR -> SE (SE of log OR, as BigaGWAS expects with OR);
             keep OR, drop LogOR and Freq
  rheuma     has no SE: SE = |ln OR| / z(P), the SE implied by OR and P.
             OR is rounded to 2 decimals, so 10% of variants read OR = 1.00
             (|ln OR| < 0.005). They are KEPT with SE = 0.005 / z(P), the SE at
             the rounding bound: dropping them shifts the median Z to -0.13 and
             LDSC munge then rejects the file. Both tools give these variants
             Z = 0 (sign of ln 1 is 0), so this SE does not affect any estimate.
             Only variants with P == 1 (z = 0, SE undefined) are dropped.

Every file is written as tab-separated, gzipped, alleles upper-cased, columns:
  SNP [CHR BP] A1 A2 (BETA|OR) SE P N

Usage:  prepare_upload.py <vignette_dir> <out_dir>
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)


def read(name):
    return pd.read_csv(src / f"{name}.sumstats.txt", sep=r"\s+", dtype={"SNP": str})


def write(df, name):
    df = df.copy()
    # upper-case alleles: the depression file uses a/c/g/t, which BigaGWAS'
    # harmonizer may fail to match against its reference
    df["A1"], df["A2"] = df["A1"].str.upper(), df["A2"].str.upper()
    path = out / f"{name}_biga.txt.gz"
    df.to_csv(path, sep="\t", index=False, float_format="%.6g",
              compression={"method": "gzip", "mtime": 0})
    print(f"{name}: {len(df)} variants, columns {list(df.columns)} -> {path.name}")


bmi = read("bmi").rename(columns={"SNPID_UKB": "SNP", "NMISS": "N"})
write(bmi[["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]], "bmi")

neuro = read("neuro").rename(columns={"POS": "BP", "N_analyzed": "N"})
write(neuro[["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]], "neuro")

dep = read("depression").rename(columns={"StdErrLogOR": "SE"})
write(dep[["SNP", "A1", "A2", "OR", "SE", "P", "N"]], "depression")

rh = read("rheuma")
z = stats.norm.isf(rh["P"] / 2)
log_or = np.abs(np.log(rh["OR"]))
log_or = np.where(log_or == 0, 0.005, log_or)        # OR printed as 1.00: rounding bound
with np.errstate(divide="ignore", invalid="ignore"):
    rh["SE"] = log_or / z
keep = np.isfinite(rh["SE"]) & (rh["SE"] > 0)
print(f"rheuma: {(rh['OR'] == 1).sum()} variants with OR == 1.00 kept (SE at rounding bound); "
      f"dropped {(~keep).sum()} with P == 1 (SE undefined)")
write(rh.loc[keep, ["SNP", "A1", "A2", "OR", "SE", "P", "N"]], "rheuma")
