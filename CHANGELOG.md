# mauricioroza/gencorr: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - unreleased

Initial version.

### Added

- Harmonization of GWAS summary statistics from common layouts into one format, with per-trait QC reports.
- LDSC 1.0.1: munge, observed + liability h2, batched cross-trait rg, sample-overlap matrix.
- LAVA 0.1.5 (`e729a24`): univariate + bivariate local rg, loci split into chunks, per-(locus, pair) seeds, `pairwise`/`joint` SNP sets, LDSC-based sample-overlap correction.
- One standard results table for all methods, validated on merge; explicit `status` for missing estimates.
- Pair modes `all`, `target`, `pairs`.
- `--min_info` / `--min_maf`: optional variant QC in harmonization, applied to all methods (LDSC munge rule: remove value <= threshold; missing kept; warning when the column is absent).
- `--lava_snp_list`: restrict LAVA to a SNP-ID list (e.g. `w_hm3.snplist`), as BigaGWAS does; recorded in `extra.snp_set`.
- LAVA: loci without shared SNPs are skipped up front (`no_snps_in_locus`), and each analysed locus logs a progress line with timing; results unchanged.
- Per-trait ancestry: samplesheet column `ancestry` (`EUR` default, `EAS`) with `--ldsc_ld_eas`, `--ldsc_w_ld_eas`, `--lava_ref_eas`. Same-ancestry pairs run LDSC/LAVA with their ancestry's references (joint LAVA groups per ancestry); cross-ancestry pairs are kept with status `cross_ancestry`. References are required only for ancestries that have same-ancestry pairs.
- rsID mapping for summary statistics without rsIDs (e.g. GRCh38 `chr:pos:ref:alt` IDs, or no SNP column): samplesheet column `build` (`GRCh37`/`GRCh38`) and `--rsid_map_grch37` / `--rsid_map_grch38` (PLINK `.bim` or CHR/BP/SNP/A1/A2 table); matched by chromosome, position and alleles, no liftover needed. Without a map, the QC report warns when most IDs are not rsIDs.
- `bin/compare_results.py` for comparison with BigaGWAS or other runs.
- BigaGWAS comparison kit (`tests/bigagwas/`, `docs/bigagwas_comparison.md`); LDSC and LAVA reproduce BigaGWAS jobs 966, 967, 971 and 972.
- Tests: Python unit + oracle tests against LDSC and LAVA; nf-test module tests for all 8 local modules (real-data + stub), function and pipeline tests.

### Changed

- Local modules follow the nf-core module specification: `task.ext.args` in every module, a `versions_<tool>` topic output for every tool used, `meta.yml` in the current nf-core format, test tags for modules run in `setup`. `nf-core modules lint` and `nf-core pipelines lint` report no failures.
- The pipeline test snapshots `pipeline_info/gencorr_software_versions.yml`.
- Module conda environments pin `python=3.12.3`, the Python the `gencorr-tools` image ships (was 3.12.8), so conda and container runs report the same versions.

### Fixed

- `tests/data/` (test LD scores, samplesheet, loci) was never committed because `.gitignore` ignored every `data/` folder; `-profile test` and the nf-tests failed on a fresh clone. The rule now ignores only a top-level `/data/`.

- `gencorr.ldsc_h2.tsv`: `samp_prev` / `pop_prev` of binary traits were NA when a trait's `.h2_liab.log` happened to be read before its `.h2.log` (file order is not guaranteed).
