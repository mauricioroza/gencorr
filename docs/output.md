# Outputs

```
results/
├── gencorr_results.tsv            all methods, one standard table (below)
├── harmonized/
│   ├── <id>.qc.json               harmonization report per trait
│   └── <id>.harmonized.tsv.gz     only with --save_harmonized
├── ldsc/
│   ├── gencorr.ldsc.tsv           LDSC rows of the standard table
│   ├── gencorr.ldsc_h2.tsv        per-trait h2 (observed + liability), intercept, ratio, lambda GC
│   ├── gencorr.sample_overlap.txt sampling-correlation matrix given to LAVA
│   ├── logs/                      raw LDSC logs: <id>.h2.log, <id>.h2_liab.log, <trait1>__bN.rg.log
│   └── munged/                    munge logs (+ .sumstats.gz with --save_intermediates)
├── lava/
│   ├── gencorr.lava.tsv           LAVA rows of the standard table
│   ├── gencorr.lava_univ.tsv      local univariate h2 and p per locus and trait
│   └── chunks/                    per-chunk tables, only with --save_intermediates
└── pipeline_info/
    ├── gencorr_software_versions.yml   exact tool versions used by every process
    └── execution_{report,timeline,trace}_*.html/txt, params_*.json
```

## The standard results table

Defined once in `bin/gencorr_common.py` (`RESULT_COLUMNS`) and enforced by
`merge_results.py`. A table that does not follow it stops the pipeline.

| column | LDSC | LAVA |
|---|---|---|
| `method` | `ldsc` | `lava` |
| `scope` | `global` | `local` |
| `trait1`, `trait2` | IDs from the samplesheet (pairs are unordered; rg is symmetric) | same |
| `locus`, `chr`, `start`, `stop` | NA | locus ID and coordinates from the locus file |
| `n_snps` | SNPs in the rg regression ("SNPs with valid alleles") | SNPs in the locus after alignment |
| `rg` | genetic correlation | local rho |
| `rg_se` | LDSC SE | NA (LAVA reports CIs, not SEs) |
| `rg_ci_low`, `rg_ci_high` | rg ± 1.96·SE (computed here) | LAVA 95% simulation CI |
| `z` | LDSC Z | NA |
| `p` | LDSC p | LAVA simulation p (bivariate) |
| `h2_1`, `h2_2` | observed-scale SNP-h2 of each trait from its own `--h2` run | local observed h2 (`run.univ`) |
| `h2_1_se`, `h2_2_se` | LDSC SE | NA |
| `status` | see below | see below |
| `extra` | JSON: `gcov_int`, `gcov_int_se`, pair-specific h2/intercept of trait2, LDSC error text, log name | JSON: `univ_p_1/2`, `h2_latent_1/2`, `r2` and its CI, `K` (#PCs), `univ_thresh`, `seed`, `snp_set` (`all` or the `--lava_snp_list` file name) |

Read it in R with `read.delim("gencorr_results.tsv", quote = "")` or in Python
with `pandas.read_csv(..., sep="\t", quoting=3)`. Parse `extra` with
`jsonlite::fromJSON` / `json.loads`.

### `status` values

| status | method | meaning |
|---|---|---|
| `ok` | both | estimate available |
| `rg_not_estimable` | ldsc | LDSC returned NA (usually h2 ≤ 0 for one trait) |
| `ldsc_error` | ldsc | LDSC raised an error for this pair; message in `extra.error` (e.g. incompatible alleles without `--ldsc_merge_alleles`) |
| `univ_filtered` | lava | one or both traits have local univariate p ≥ `--lava_univ_thresh`; not tested (LAVA's standard practice) |
| `estimate_out_of_bounds` | lava | LAVA set the estimate to NA because it exceeded `--lava_param_lim` |
| `trait_dropped_in_locus` | lava | LAVA dropped a trait in this locus (invalid local estimates) |
| `no_snps_in_locus` | lava | no SNP shared by the traits and the reference falls in the locus; skipped without calling LAVA (fast) |
| `locus_not_processed` | lava | LAVA could not process the locus (e.g. too few SNPs or PCs) |
| `bivar_error` | lava | `run.bivar` raised an error (message in the task log) |
| `cross_ancestry` | both | the two traits have different `ancestry` in the samplesheet; LDSC and LAVA assume one population, so the pair is not estimated. One row per method (locus NA); `extra` holds `ancestry_1` / `ancestry_2` |

## `lava/gencorr.lava_univ.tsv`

One row per (locus, trait) that LAVA processed:

| column | meaning |
|---|---|
| `locus`, `chr`, `start`, `stop` | locus ID and coordinates |
| `phenotype` | trait ID |
| `h2_obs` | local observed-scale h2 (`run.univ`) |
| `h2_latent` | local liability-scale h2; only for binary traits with `pop_prev` (else NA) |
| `p` | univariate p-value (compared with `--lava_univ_thresh`) |
| `n_snps`, `K` | SNPs in the locus and number of PCs LAVA kept |
| `loaded_traits` | traits loaded together in that LAVA call (`;`-separated): the pair with `pairwise`, all traits with `joint` |

## Precision

* LDSC values are copied from LDSC's logs, which print 4 decimals. That is
  the precision of all LDSC numbers here and in BigaGWAS.
* LAVA values are written with R's default 15 significant digits. CIs and
  p-values are Monte Carlo estimates (1e4 to 1e6 draws). They are exactly
  reproducible for a given `--lava_seed`, but a different seed gives slightly
  different values.

## `harmonized/<id>.qc.json`

`rows_in`, `rows_out`, the detected `column_mapping`, `z_source_counts` (how
many Z-scores came from a Z column / P / BETA-SE), `removed` (per filter),
`p_recomputed_from_z`, `lambda_gc`, `mean_chi2`, `median_n`, the QC
`thresholds` used (`min_info`, `min_maf`), where MAF came from (`maf_source`),
how many variants without INFO/MAF were kept (`info_missing_kept`,
`maf_missing_kept`), and warnings such as `warning_effective_n`,
`warning_min_info`, `warning_min_maf`, `warning_non_rsid` (most SNP IDs are
not rsIDs and no rsID map was given). Variants removed by the thresholds are
counted in `removed` as `info_le_<t>` / `maf_le_<t>`. With an rsID map,
`rsid_mapping` gives the map file, the build, `rows_without_rsid`, `mapped`,
`unmatched` and how variants were matched, and `removed.no_rsid_after_mapping`
counts the variants dropped because no rsID was found. Look at these first when a trait behaves
unexpectedly.
