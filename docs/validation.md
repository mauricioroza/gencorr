# Validation and tests

The question the tests answer is *"does the pipeline compute what LDSC and LAVA
compute?"*. The trusted reference is always the original tool, run directly.
You don't need to read the pipeline code to trust it: the comparisons are
short, readable scripts.

## Test layers

| layer | command | what it proves | runtime |
|---|---|---|---|
| unit (Python) | `pytest tests/python -k "not oracle"` | harmonizer rules (column detection, Z derivation, filters, N fallback, byte-identical reruns); LDSC log parsing matches the numbers in real logs; schema checks; comparison tool | 20 s |
| **LDSC oracle** | `pytest tests/python/test_ldsc_oracle.py` | raw file → `munge_sumstats.py` → `ldsc --rg` (what you or BigaGWAS do by hand) gives **identical** rg/se/z/p/intercepts to the pipeline route; batching pairs into one `--rg` call changes nothing | 2 min |
| **LAVA oracle** | `pytest tests/python/test_lava_oracle.py` | `bin/run_lava.R` gives **identical** rho, CIs, p and local h2 to LAVA called as in its vignette (`tests/oracle/lava_native.R`); splitting loci into chunks changes nothing; harmonized input vs LAVA reading raw files: same SNPs/statuses, rg within 1e-4 (explained below) | 15 min |
| function (nf-test) | `nf-test test subworkflows/local/utils_nfcore_gencorr_pipeline/tests/` | pair construction for all/target/pairs modes, error cases | 40 s |
| pipeline (nf-test) | `nf-test test tests/default.nf.test --profile +apptainer` | full test profile succeeds; row counts; snapshot of every (method, pair, locus, status, rg to 3 decimals); `-stub` wiring; target + joint mode | 15 min |

**Nextflow ≥ 26 and nf-test:** the pipeline itself is written in Nextflow's
strict syntax (`nextflow lint .` reports no errors) and runs on 26.x as is. nf-test
0.9.3, however, generates wrapper scripts with syntax the strict parser rejects.
Run nf-test with the legacy parser:

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/default.nf.test --profile +docker
```

Oracle tests need `ldsc.py`/`munge_sumstats.py` and R with LAVA on `PATH`, and
are skipped otherwise. In the container they are available:

```bash
apptainer exec gencorr-tools_0.1.0.sif pytest tests/python   # LAVA + unit tests
apptainer exec docker://quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2 ...  # LDSC
```

## Status of the test suite at release (v0.1.0dev)

Run on Linux, with Nextflow 25.10.6, nf-test 0.9.3, LDSC 1.0.1 (bioconda,
Python 2.7.15), LAVA 0.1.5 (`e729a24`) on R 4.3.3:

* Python: **41/41 pass** (unit + LDSC oracle + LAVA oracle)
* nf-test function tests: **6/6 pass**
* nf-test pipeline tests: **3/3 pass** on Nextflow 25.10 and 26.10.5 (`-profile test`, `-stub`, target + joint);
  snapshot in `tests/default.nf.test.snap`
* `-profile test` also run on Nextflow 26.10.5 (strict syntax): results identical
  to the 25.10 run at zero tolerance (66/66 rows)
* Tools were on PATH (`-profile test,local_tools`); nf-schema plugin 2.5.1 was
  used locally because the plugin registry was unreachable from the build machine

Not yet run: the pipeline inside the published containers (`-profile
apptainer/docker`), because the images had not been built at the time. **Run
`-profile test,apptainer` once after building the image**; it must produce
the same `tests/default.nf.test.snap`.

## Known, explained differences

| comparison | difference | cause |
|---|---|---|
| LAVA on harmonized vs raw files | rg differs in the 5th decimal, local h2 by ~1e-4 relative, for traits given only as OR/LogOR | variants with OR exactly 1 (435 in the vignette depression file): LAVA's reader gives them sign −1 (`ifelse(OR > 1, 1, -1)`); the harmonizer gives Z = 0 because there is no direction |
| LDSC global rg vs BigaGWAS | expected identical when the same munged inputs, LD scores and `w_hm3` list are used | same LDSC version and options |
| LAVA vs BigaGWAS | depends on BigaGWAS settings that are not in its public repo (`specific_lava.R`, `2traits.R`, `lava_pgc.R`) | see below |

## Comparing with BigaGWAS

`bin/compare_results.py` matches pairs order-insensitively, applies
per-column tolerances, and exits non-zero on failure (usable in CI).

### LDSC

1. Download the CSV of a BigaGWAS LDSC job (`Trait1, Trait2, rg, se, z, p_value, ...`).
2. Run gencorr on the same traits with `--pair_mode target --targets <trait1>`
   and the same LD scores (`--ldsc_ld`, `--ldsc_merge_alleles w_hm3.snplist`).
3. Map BigaGWAS names to your IDs (`biga_names.csv`: `reference_name,gencorr_id`) and compare:

```bash
bin/compare_results.py --test results/gencorr_results.tsv --method ldsc \
    --reference biga_ldsc.csv --reference-format bigagwas_ldsc \
    --trait-map biga_names.csv --report ldsc_vs_biga
```

Default tolerances: rg/se 2e-4 (two 4-decimal roundings), p 0.01 on the log10 scale.

### LAVA

Set the same parameters BigaGWAS used: `--lava_univ_thresh`,
`--lava_param_lim`, CIs/p on, the same LD reference and loci, and **ask the
BigaGWAS team** which of these they used (not in their public code):

* `--lava_snp_set joint` (all group traits loaded together) or `pairwise`
* whether a sample-overlap matrix was used (BigaGWAS has an `overlap` option for 2-trait jobs)
* the LAVA version installed on their server

```bash
bin/compare_results.py --test results/gencorr_results.tsv --method lava --only-ok \
    --reference biga_lava.tsv --reference-format custom \
    --colmap trait1=phen1,trait2=phen2,locus=locus,rg=rho,p=p \
    --tol rg=1e-3,p=0.2 --allow-missing --report lava_vs_biga
```

With a different seed, LAVA's Monte Carlo p-values and CIs differ slightly
between any two runs, so compare p on the log scale with a loose tolerance.
`rg` (rho) is not simulated and should agree closely when settings match.

## Regenerating test data

`tests/data/ldsc/` (LD scores + allele list for the test reference) is built by
`tests/data/make_test_data.sh` (needs ldsc and plink2). Everything else in the
test profile is downloaded from the pinned LAVA commit at run time.
