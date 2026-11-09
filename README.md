# mauricioroza/gencorr

Global and local genetic correlation from GWAS summary statistics, as a
Nextflow pipeline built on the nf-core template.

* **Global rg**: LDSC 1.0.1 (Bulik-Sullivan et al. 2015), same version BigaGWAS uses
* **Local rg**: LAVA 0.1.5 (Werme et al. 2022), pinned to commit `e729a24`
* **Pluggable**: every method reads the same harmonized inputs and writes the
  same results table, so adding HDL, SumHer, Popcorn, HDL-L … is one subworkflow
  ([docs/extending.md](docs/extending.md)).
* **Testable**: 43 Python tests + 22 nf-test module tests + 6 function tests + 3 pipeline tests,
  including *oracle* tests that check the pipeline reproduces LDSC and LAVA
  run directly ([docs/validation.md](docs/validation.md)).
* **Checked against BigaGWAS**: on the same data and settings, LDSC rg is
  identical and LAVA rho agrees to 1e-5–1e-7
  ([docs/bigagwas_comparison.md](docs/bigagwas_comparison.md#6-results-on-the-test-jobs)).

```
 samplesheet.csv ─► HARMONIZE_SUMSTATS ─┬─► METHOD_LDSC ── munge ─ h2 ─ rg ─► ldsc table ─┐
 (any common GWAS      (one format,      │                          └─► sample overlap ──┐ │
  column layout)        QC report)       │                                               ▼ │
                                         └─► METHOD_LAVA ── loci split into chunks ─► lava table
                                                                                          │
                                         gencorr_results.tsv  ◄── MERGE (schema-checked) ◄┘
```

## Quick start

```bash
# 1. Check your setup on the bundled test data (4 traits, 6 pairs, ~5 min)
nextflow run mauricioroza/gencorr -profile test,apptainer --outdir test_out

# 2. Real run
nextflow run mauricioroza/gencorr -profile apptainer \
    --input samplesheet.csv --outdir results \
    --methods ldsc,lava \
    --ldsc_ld /refs/eur_w_ld_chr --ldsc_merge_alleles /refs/w_hm3.snplist \
    --lava_ref /refs/lava-ukb/lava-ukb-v1.1
```

`samplesheet.csv` has one row per trait:

```csv
id,sumstats,n,n_cases,n_controls,pop_prev,colmap,label
bmi,/data/bmi.txt.gz,,,,,,Body mass index
mdd,/data/mdd.tsv.gz,,170756,329443,0.15,,Major depression
scz,/data/scz.txt,,53386,77258,0.01,"Z=STAT",Schizophrenia
```

Column names in the files are detected automatically (SNP/rsid/MarkerName,
A1/effect_allele, BETA/OR/LogOR, SE, P, Z, N/NMISS/…); `colmap` overrides
when needed. See [docs/usage.md](docs/usage.md) for pair modes (`all`,
`target` = BigaGWAS "one vs group", `pairs`) and reference-data downloads.

## Outputs

`results/gencorr_results.tsv`: one row per (method, trait pair, locus); main
columns shown (all columns in [docs/output.md](docs/output.md)):

| method | scope | trait1 | trait2 | locus | rg | rg_se | rg_ci_low | rg_ci_high | p | h2_1 | h2_2 | status | extra |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ldsc | global | bmi | mdd | NA | 0.11 | 0.07 | … | … | 0.11 | … | … | ok | {"gcov_int": …} |
| lava | local | bmi | mdd | 464 | 0.51 | NA | 0.27 | 0.81 | 1.1e-4 | … | … | ok | {"r2": …} |
| lava | local | bmi | mdd | 100 | NA | NA | NA | NA | NA | … | … | univ_filtered | … |

Pairs that could not be estimated are **kept with a status** (`univ_filtered`,
`ldsc_error`, `estimate_out_of_bounds`, …) instead of disappearing.
Details: [docs/output.md](docs/output.md).

## Documentation

| | |
|---|---|
| [docs/usage.md](docs/usage.md) | samplesheet, parameters, reference data, running on UPPMAX |
| [docs/output.md](docs/output.md) | every output file and column, status values |
| [docs/validation.md](docs/validation.md) | the test layers, how to compare against BigaGWAS |
| [docs/bigagwas_comparison.md](docs/bigagwas_comparison.md) | running the same data on BigaGWAS and comparing |
| [docs/extending.md](docs/extending.md) | adding a new method, step by step |
| [docs/design.md](docs/design.md) | design decisions and where results differ from BigaGWAS on purpose |

## Requirements

* Nextflow ≥ 25.10.4
* Apptainer/Singularity or Docker. LAVA has no conda package, so `-profile conda`
  works for LDSC only (`--methods ldsc`).
* The helper image `ghcr.io/mauricioroza/gencorr-tools:0.1.0` must be built once
  from [containers/gencorr-tools](containers/gencorr-tools/Dockerfile)
  (instructions in the Dockerfile). Note: LAVA's licence says *All rights
  reserved*, so keep that image private unless the LAVA authors agree.

## Credits

Written by MR (Genis). Method choices
follow BigaGWAS (bigagwas.org), whose LDSC/LAVA usage this pipeline reproduces
and extends. Tool citations: [CITATIONS.md](CITATIONS.md).
