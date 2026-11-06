# mauricioroza/gencorr: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - unreleased

Initial version.

### Added

- `--min_info` / `--min_maf`: optional variant QC in harmonization, applied to all methods (LDSC munge rule: remove value <= threshold; missing kept; warning when the column is absent).
- `--lava_snp_list`: restrict LAVA to a SNP-ID list (e.g. `w_hm3.snplist`), as BigaGWAS does; recorded in `extra.snp_set`.
- BigaGWAS comparison kit (`tests/bigagwas/`, `docs/bigagwas_comparison.md`); LDSC and LAVA reproduce BigaGWAS on the test data.

- Harmonization of GWAS summary statistics from common layouts into one format, with per-trait QC reports.
- LDSC 1.0.1: munge, observed + liability h2, batched cross-trait rg, sample-overlap matrix.
- LAVA 0.1.5 (`e729a24`): univariate + bivariate local rg, loci split into chunks, per-(locus, pair) seeds, `pairwise`/`joint` SNP sets, LDSC-based sample-overlap correction.
- One standard results table for all methods, validated on merge; explicit `status` for missing estimates.
- Pair modes `all`, `target`, `pairs`.
- `bin/compare_results.py` for comparison with BigaGWAS or other runs.
- Tests: Python unit + oracle tests against LDSC and LAVA, nf-test function and pipeline tests.
