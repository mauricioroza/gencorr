# Design decisions

Why the pipeline is built the way it is, and where it deliberately departs
from the tools' default behaviour.

## Why Nextflow + existing tools, not a Rust rewrite

The statistics happen inside LDSC (Python 2), LAVA (R + C++), LDAK and
Popcorn, not in the orchestration. Rewriting the orchestration in a faster
language would not speed up the computation. Rewriting LDSC/LAVA themselves means re-validating
statistical methods, and LAVA's licence (*All rights reserved*) does not permit
derivative implementations without permission.

The speed comes from parallelism instead:

* LAVA loci are independent, so they are split into chunks and run as separate
  jobs (`--lava_loci_per_chunk`). Wall time scales with the number of nodes.
* LDSC munging and h2 run once per trait and are cached (`-resume`); rg runs in
  batches per trait1.

If profiling (`pipeline_info/execution_trace_*.txt`) ever shows one step
dominating, that single step can be replaced by a compiled tool (e.g. the Rust
LDSC port, `crates.io/crates/ldsc`) behind the same module interface, and the
oracle tests will show whether it reproduces the original.

## One harmonization step for all methods

One script produces one harmonized file per trait,
with a QC report, and each method only *reformats* it
(`export_sumstats.py`). All methods therefore see the same variants, alleles
and Z-scores, and differences between methods are methodological, not
caused by input handling.

## Departures from the tools' defaults (all documented in tests)

| topic | here | why |
|---|---|---|
| LAVA SNP intersection | `pairwise` by default: each pair uses SNPs shared by its two traits | LAVA intersects SNPs across *all* loaded traits. A pair's result then depends on which other traits happen to be in the job. `--lava_snp_set joint` restores that behaviour. |
| LAVA empty loci | skipped before `process.locus()` (`no_snps_in_locus`) | LAVA scans the full reference SNP table for every locus (~1 s each with 1000G EUR) even when the locus has no data; on sparse input this was ~40 min of a ~45 min run. Results for all other loci are identical (tested) |
| LAVA bivariate call | `run.bivar` on exactly two traits | with more traits LAVA draws the CIs for all pairs from one joint simulation, which ties a pair's CI to the other traits |
| LAVA seeds | seed per (locus, trait1, trait2) | results identical regardless of chunking, job order or other pairs |
| LAVA sample overlap | from LDSC intercepts by default (LAVA vignette recipe, symmetric average of both intercept orders) | sample overlap biases local rg upwards |
| pairs without estimate | kept with a `status` | nothing silently disappears; the status says why there is no estimate |
| LDSC liability h2 | second `--h2` run with prevalences; both scales reported from LDSC itself | LDSC prints only the liability value when prevalences are given |
| OR = 1 variants | Z = 0 | LAVA's raw reader gives them sign −1 |
| gzip outputs | written with mtime 0 | byte-identical reruns, stable checksums |
| LAVA SNP list | all SNPs by default; HapMap3 (or any list) with `--lava_snp_list` | the default follows the LAVA authors, who do not filter; a HapMap3 list keeps fewer SNPs (~16k of ~100k on the test data) |
| LAVA version | 0.1.5 | current release at the time of writing; 0.1.5 changed the model for binary traits |
| LAVA case/control | used when `n_cases`/`n_controls` are given | binary traits use LAVA's case/control model instead of being treated as continuous |
| SumHer version | LDAK 6.1 | LDAK 6 writes `SE` where 5.x wrote `SD` (the parser accepts both); its default without `--genomic-control` models multiplicative inflation, so gencorr passes `--genomic-control NO` explicitly (`--sumher_genomic_control`) |
| SumHer Z | the harmonized Z (a Z column, else from P, by default; `--z_priority`) for every method | all methods see the same Z; `--z_priority beta_se,p` gives Z = BETA/SE instead |
| Popcorn p-value | `p` tests rg ≠ 0 (z = rg/SE), like the other methods; Popcorn's own rg < 1 test kept in `extra` | Popcorn reports the p for rg < 1; mixing the two in one column would make cross-ancestry results look significant for the opposite reason |
| Popcorn input | harmonized Z only (no beta/SE), so every method uses the same Z | Popcorn recomputes Z from beta/SE when both are given; on Popcorn's test data the two routes agree to 9 digits |
| ancestry | per trait (`ancestry` column); each same-ancestry pair uses its ancestry's references; cross-ancestry pairs kept as `cross_ancestry` | mixed-ancestry sets run in one job; LAVA runs for EAS too when an EAS reference is given |

## Reproducibility

* Every tool version is pinned: LDSC by bioconda build, LAVA by Git commit,
  R packages by a dated CRAN snapshot (rocker), Python packages by `==`.
  Python itself is the image's Ubuntu 24.04 `python3` (3.12.3), and the
  modules' `environment.yml` files pin the same versions, so conda and
  container runs report identical software versions.
* `pipeline_info/gencorr_software_versions.yml` records the versions actually
  used by each process. `params_*.json` records the parameters.
* Test data are fetched from a pinned LAVA commit; derived LD scores are
  regenerated by a script.

## Things intentionally not done (yet)

* **Liftover**: not needed. Every method matches variants by rsID and takes
  positions from its own reference, so files with rsIDs work on any build.
  Files without rsIDs get them from an rsID map on their own build
  (`--rsid_map_grch37/38`), by chromosome, position and alleles. Lifting
  every file between builds and assigning rsIDs from dbSNP would be slow and
  is not needed.
* **Allele alignment to a reference panel** in the harmonizer: LDSC aligns via
  `--merge-alleles`, and LAVA aligns to its LD reference itself.
* **Multivariate LAVA** (partial correlation, multiple regression): the runner
  is structured so `run.pcor`/`run.multireg` could be added as another output
  type once there is a use case.
