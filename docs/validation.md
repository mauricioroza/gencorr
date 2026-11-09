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
| module (nf-test) | `nf-test test --tag modules_local` | each of the 8 local modules on its own, with small real inputs: output layout, the numbers the tool reports (LDSC summary lines; LAVA status, n SNPs, K and rho per locus), error paths (duplicated estimates rejected), and the `stub:` block | 10 min |
| function (nf-test) | `nf-test test subworkflows/local/utils_nfcore_gencorr_pipeline/tests/` | pair construction for all/target/pairs modes, error cases | 40 s |
| pipeline (nf-test) | `nf-test test tests/default.nf.test --profile +apptainer` | full test profile succeeds; row counts; snapshot of every (method, pair, locus, status, rg to 3 decimals) and of the software versions; `-stub` wiring; mixed-ancestry routing (same-ancestry pairs estimated, cross-ancestry pairs `cross_ancestry`, joint LAVA groups per ancestry); target + joint mode | 20 min |
| nf-core lint | `nf-core pipelines lint` and `nf-core modules lint --all --local` | template conformance; every module's `meta.yml` matches its channels (official schema) | 1 min |

### Running the Python tests

The scripts in `bin/` are started as `python3 …`, so the environment with
pandas/numpy/scipy must be **first on `PATH`**, not only the one running pytest.
Use a virtual environment with the container's pins:

```bash
python3 -m venv ~/venvs/gencorr
~/venvs/gencorr/bin/pip install pandas==2.2.3 numpy==2.2.2 scipy==1.15.1 pytest nf-core
source ~/venvs/gencorr/bin/activate      # every new shell
pytest tests/python                      # oracle tests skip unless LDSC / LAVA are on PATH
```

**Module tests** live next to each module (`modules/local/<name>/tests/main.nf.test`
and `.snap`), as in nf-core/modules. Their inputs are the LAVA vignette files
(pinned commit, `params.lava_test_data` in `tests/nextflow.config`),
`tests/data/` and the LDSC log fixtures in `tests/python/fixtures/ldsc`.
Modules that need upstream output build it in a `setup` block (for example
`LDSC_RG` runs harmonize → export → munge first). Run a single module with
its tag:

```bash
NXF_SYNTAX_PARSER=v1 nf-test test --tag ldsc/rg --profile +docker
```

Snapshots hold only values that are deterministic across machines: md5 of our
own outputs (nf-test hashes `.gz` files on their uncompressed content), LDSC's
printed summary lines, and LAVA rho rounded to 3 decimals. Module snapshots
record only the tool's own version (`ldsc 1.0.1`, `lava 0.1.5`), so they hold
with conda, Docker or tools on `PATH`. The **pipeline** snapshot also records
`pipeline_info/gencorr_software_versions.yml` (nf-core template), which includes
the Python and R builds; it holds for Docker/Apptainer and for conda because the
`environment.yml` files pin exactly what the image ships (Python 3.12.3,
pandas 2.2.3, numpy 2.2.2, scipy 1.15.1, R 4.4.2), but not with arbitrary tools
on `PATH`. The LAVA module test runs without sample overlap, so its rho values
differ from the pipeline snapshot, which uses the LDSC overlap matrix.

**Nextflow ≥ 26 and nf-test:** the pipeline itself is written in Nextflow's
strict syntax (`nextflow lint .` reports no errors) and runs on 26.x as is. nf-test
(checked with 0.9.3–0.9.5), however, generates wrapper scripts with syntax the
strict parser rejects. Run nf-test with the legacy parser:

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/default.nf.test --profile +docker
```

nf-test's own `-c` option is the *nf-test* config file, not a Nextflow config,
so it cannot pass `local.config`; the test config (`tests/nextflow.config`) is
used instead.

**Test data download:** the test profile fetches the LAVA vignette files from
`raw.githubusercontent.com` on every run. When GitHub throttles anonymous
downloads, Nextflow waits indefinitely while staging them and the test appears
to hang. Check the speed first, e.g.
`curl -o /dev/null -w '%{speed_download}\n' <one of the URLs in tests/data/samplesheet.csv>`.

Oracle tests need `ldsc.py`/`munge_sumstats.py` and R with LAVA on `PATH`, and
are skipped otherwise. The tools are in the containers, but **the
`gencorr-tools` image does not contain pytest** (and the LDSC biocontainer is
Python 2), so the oracle tests cannot yet be run inside them as is. Options:
add `pytest` to `containers/gencorr-tools/Dockerfile` at the next image
version, or run them on a machine with LDSC and LAVA installed.

## Status of the test suite (v0.1.0dev)

**2026-11-09, in Docker** (`--profile +docker`, `gencorr-tools:0.1.0` + LDSC
biocontainer), Nextflow 26.10.5, nf-test 0.9.5, WSL2 Ubuntu:

* nf-test module tests: **23/23 pass** (incl. rsID mapping of a GRCh38-style
  file); function tests **6/6**; pipeline tests **4/4** (`-profile test`,
  `-stub`, mixed EUR + EAS ancestry, target + joint), all with `--ci`
* Python: **48 pass, 6 skipped** (oracle tests: LDSC/LAVA not on `PATH`)
* `nf-core pipelines lint`: 0 failures; `nf-core modules lint --all --local`: 0 failures
  (remaining warnings are expected for a private pipeline: private container
  registry, custom meta keys, underscores in local module names)

**Earlier (first build), tools on `PATH`** (`-profile test,local_tools`),
Nextflow 25.10.6, nf-test 0.9.3, LDSC 1.0.1 (bioconda, Python 2.7.15),
LAVA 0.1.5 (`e729a24`) on R 4.3.3:

* Python: **41/41 pass** including both oracle layers (the 2 tests added since
  are unit tests)
* nf-test module, function and pipeline tests pass on Nextflow 25.10 and 26.10.5;
  module snapshots identical with pandas 2.2.3 and 3.0.5
* `-profile test` on 25.10 and 26.10.5: identical results at zero tolerance (66/66 rows)

Not yet run: Apptainer/Singularity (`-profile test,apptainer`) and the `conda`
profile (LDSC only), and nf-test in CI with the private image.

## Known, explained differences

| comparison | difference | cause |
|---|---|---|
| LAVA on harmonized vs raw files | rg differs in the 5th decimal, local h2 by ~1e-4 relative, for traits given only as OR/LogOR | variants with OR exactly 1 (435 in the vignette depression file): LAVA's reader gives them sign −1 (`ifelse(OR > 1, 1, -1)`); the harmonizer gives Z = 0 because there is no direction |
| LDSC global rg vs BigaGWAS | identical on job 972; 4th-decimal difference on job 966 | same LDSC version and options; BigaGWAS' own harmonization can drop a few SNPs (966: 4) |
| LAVA vs BigaGWAS | rho within 1e-7 to 1e-5 with `--lava_snp_list w_hm3.snplist`; p/CIs differ by Monte Carlo noise | see [bigagwas_comparison.md](bigagwas_comparison.md#6-results-on-the-test-jobs) |

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
`--lava_param_lim`, CIs/p on, the same LD reference and loci. From the job logs
of the test jobs (967, 971), BigaGWAS' **2-trait** LAVA form:

* restricts LAVA to a HapMap3-based SNP list → use `--lava_snp_list w_hm3.snplist`
* treats **all traits as continuous** ("All phenotypes treated as CONTINUOUS")
  → leave `n_cases`/`n_controls` empty in the samplesheet
* uses sample overlap 0 unless set in the job → `--lava_overlap none`
* runs LAVA 0.1.0 (gencorr: 0.1.5; no visible effect for continuous traits)

Still unknown, because the scripts are not in BigaGWAS' public code
(`specific_lava.R`, `2traits.R`, `lava_pgc.R`): the settings of its group
("massive") LAVA jobs, e.g. whether all group traits are loaded together
(`--lava_snp_set joint`).

```bash
bin/compare_results.py --test results/gencorr_results.tsv --method lava --only-ok \
    --reference biga_lava.tsv --reference-format custom \
    --colmap trait1=phen1,trait2=phen2,locus=locus,rg=rho,p=p \
    --tol rg=1e-3,p=0.2 --allow-missing --report lava_vs_biga
```

BigaGWAS writes `<job>_out_bivar.txt` as **comma**-separated; convert it to
tabs first (`sed 's/,/\t/g'`) and map its trait names with `--trait-map`.

With a different seed, LAVA's Monte Carlo p-values and CIs differ slightly
between any two runs, so compare p on the log scale with a loose tolerance.
`rg` (rho) is not simulated and should agree closely when settings match.

## Regenerating test data

`tests/data/ldsc/` (LD scores + allele list for the test reference) is built by
`tests/data/make_test_data.sh` (needs ldsc and plink2). Everything else in the
test profile is downloaded from the pinned LAVA commit at run time.
