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
- Popcorn cross-ancestry genetic correlation, `--methods popcorn`: runs on the EUR × EAS pairs (the other methods take same-ancestry pairs; Popcorn reports those as `same_ancestry`). Modules `POPCORN_FIT` and `POPCORN_SUMMARISE` (`bin/parse_popcorn.py`), subworkflow `METHOD_POPCORN`, `popcorn` export format; `--popcorn_cfile`, `--popcorn_cfile_ancestries`, `--popcorn_gen_effect`, `--popcorn_use_mle`; liability-scale h2 for binary traits. The `p` column tests rg ≠ 0; Popcorn's own rg < 1 test is kept in `extra`. Image `containers/gencorr-popcorn` (Python 3.11, pinned deps incl. numpy 1.23.5 and six, Popcorn commit 20a29e8; GPL-3.0). Failed runs kept as `popcorn_error`.
- SumHer (LDAK 6.1) global genetic correlation, `--methods sumher`: modules `SUMHER_RG` (`ldak6 --sum-cors`, one run per same-ancestry pair) and `SUMHER_SUMMARISE` (`bin/parse_sumher.py`), subworkflow `METHOD_SUMHER`, `sumher` export format. Tagging file per ancestry (`--sumher_tagfile`, `--sumher_tagfile_eas`), `--sumher_snp_positions` to name variants `CHR:BP` for tagging files with `CHR:BP` names (older LDAK releases; current pre-computed files use rsIDs), BigaGWAS' options (`--sumher_cutoff 0.01`, `--sumher_genomic_control false`). LDAK image `containers/gencorr-ldak` (GPL-3.0, publishable). Failed LDAK runs are kept as `sumher_error`. Included in `-profile test` with a small tagging file (`tests/data/sumher`).
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

- SumHer: when no SNP of a trait is in the tagging file, `extra.n_snps_1/2` is 0 for that trait and empty for a trait LDAK never checked (was the number of SNPs in the tagging file). Docs: current LDAK pre-computed tagging files use rsIDs, so `--sumher_snp_positions` is only for older `CHR:BP`-named files.

- `tests/data/` (test LD scores, samplesheet, loci) was never committed because `.gitignore` ignored every `data/` folder; `-profile test` and the nf-tests failed on a fresh clone. The rule now ignores only a top-level `/data/`.

- `gencorr.ldsc_h2.tsv`: `samp_prev` / `pop_prev` of binary traits were NA when a trait's `.h2_liab.log` happened to be read before its `.h2.log` (file order is not guaranteed).
