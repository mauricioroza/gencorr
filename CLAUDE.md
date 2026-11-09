# gencorr — notes for Claude

Nextflow (nf-core template) pipeline for global and local genetic correlation (LDSC + LAVA),
built to replace BigaGWAS (bigagwas.org) for the Genis group. Owner: GitHub `mauricioroza`,
repo `mauricioroza/gencorr` (private), working branch `dev` (also `master`, `TEMPLATE`).

Goals: reproducible, testable, well documented; results checkable against the original tools;
new methods easy to add. A Rust rewrite was considered and rejected (stats run inside LDSC/LAVA;
BigaGWAS's slowness is rsID mapping, which this pipeline avoids by design).

## Environment (WSL)

- Repo lives in WSL Ubuntu: `/home/rozam/gencorr` (Windows: `\\wsl.localhost\Ubuntu\home\rozam\gencorr`).
- Run git, nextflow, nf-test, docker, python **inside WSL** (`wsl.exe -d Ubuntu -- bash script.sh`).
  Windows-side git fails with "dubious ownership". Quoting through PowerShell → wsl breaks easily;
  write a small script and run it with bash instead.
- Editing files through the Windows `\\wsl.localhost\...` path drops the executable bit: after editing
  anything in `bin/`, run `chmod +x bin/<script>` in WSL (otherwise tasks fail with exit 126).
- nf-test's `-c` is the nf-test config, not a Nextflow config; it cannot take `local.config`.
- WSL has ~7.5 GB RAM: always run with `-c local.config` (resourceLimits cpus 4, memory 6.GB, time 4.h).
- nf-test on Nextflow ≥ 26 needs `NXF_SYNTAX_PARSER=v1` (nf-test wrappers only; pipeline doesn't), e.g.
  `NXF_SYNTAX_PARSER=v1 nf-test test --tag modules_local --profile +docker`.
- Keep code Nextflow strict-syntax clean (`nextflow lint` = 0 errors; runs on 25.10 and 26.10.5):
  no implicit `it`, no `switch`/`for` in workflow code, unused closure params prefixed `_`.

## Architecture (see docs/design.md, docs/extending.md)

- Harmonized input contract `HARMONIZED_COLUMNS`: SNP CHR BP A1 A2 Z BETA SE P N EAF INFO. Methods read only this.
- Results table `RESULT_COLUMNS`: one row per (method, trait1, trait2, locus), explicit `status`,
  method-specific values in JSON `extra` column. `MERGE_RESULTS` validates schema, rejects duplicates.
- `methodRegistry()` records each method's required references and conda support.
- Each method = one `METHOD_<NAME>` subworkflow (`take: ch_harmonized, ch_pairs, refs` / `emit: results`).
  Candidates (HDL, SumHer, Popcorn, HDL-L) and the merge checklist are in docs/extending.md.
- `bin/harmonize_sumstats.py`: column synonyms (STAT never auto-mapped), Z from z → p → beta/se,
  OR → ln OR, N from column → `--n` → cases+controls; drops non-SNVs, A1==A2, missing Z/N, all copies of
  duplicates, optionally ambiguous; per-trait QC JSON (incl. lambda_gc). Byte-reproducible output
  (gzip mtime 0, `%.8g`) — preserve this.
- LDSC 1.0.1 (bioconda, same as BigaGWAS): munge with merge-alleles (`--signed-sumstats Z,0`), h2
  (observed + liability when n_cases/n_controls/pop_prev given), batched `--rg` (`ldsc_rg_batch_size=50`),
  summary step also builds the sample-overlap matrix from intercepts (used by LAVA).
- LAVA 0.1.5 pinned at commit e729a24 (BigaGWAS uses 0.1.0 — known difference). `bin/run_lava.R`:
  loci chunked (`lava_loci_per_chunk=250`); per-(locus,trait1,trait2) seeds so results are chunk-invariant;
  adaptive Monte Carlo (1e4→1e5→1e6); `--lava_snp_set pairwise|joint`; `--lava_overlap ldsc`;
  `--lava_snp_list` (e.g. HapMap3 w_hm3.snplist, matches BigaGWAS); loci with no shared SNPs skipped as
  `no_snps_in_locus`. Statuses: ok, univ_filtered, trait_dropped_in_locus, estimate_out_of_bounds,
  bivar_error, locus_not_processed, no_snps_in_locus.
- QC: `--min_info` / `--min_maf` (0 = off; drop when value ≤ threshold, keep missing, warn if column absent).
- Pair modes: all, target, pairs (self-pairs rejected, reversed duplicates removed).

## Testing / validation (docs/validation.md, docs/bigagwas_comparison.md)

- Principle: oracle tests — the reference is always the original tool run directly.
- Python unit + oracle: 54 tests (`tests/python`; 6 oracle tests skip unless LDSC/LAVA are on PATH).
  Run with the venv active: `source ~/venvs/gencorr/bin/activate && pytest tests/python` (scripts are
  called as `python3`, so the venv must be first on PATH). nf-test: module (22), function (6),
  pipeline (3: `-profile test`, `-stub`, target + joint).
- nf-core tools 4.2.1 is in `~/venvs/gencorr`: `nf-core pipelines lint`, `nf-core modules lint --all --local`
  (run from the repo root). Module warnings about the private container, custom meta keys, underscores in
  names and `process_single`+`process_long` are expected for this private pipeline.
- Test data is downloaded from raw.githubusercontent.com on every pipeline test; GitHub can throttle it
  (seen 2026-11-07, ~150 KB/s), which makes Nextflow hang on staging. Check speed with curl first.
- nf-test snapshots must contain only machine-independent values (identical under pandas 2.2.3 and 3.0.5).
- BigaGWAS comparison kit: `tests/bigagwas/`. Jobs 966/967/971/972 agree (differences explained:
  SNP set, LAVA version, Monte Carlo noise, OR==1 sign handling in LAVA's own reader).
- BigaGWAS example data (ieu-a-25 T2D, ieu-a-9 CHD) is Metabochip, not genome-wide: plumbing tests only,
  LDSC not valid on it.

## Rules

- Never commit user data: `biga_*` dirs, `results*`, `*.tar`, `test_out/`, `work/`, `.nextflow.log*`.
- LAVA is "All rights reserved": the image `ghcr.io/mauricioroza/gencorr-tools:0.1.0` must stay private.
- Commits authored as `Mauricio Roza <113365401+mauricioroza@users.noreply.github.com>`.

## Status (as of 2026-11-08)

- On dev and pushed: module nf-tests; parse_ldsc log-order fix (samp_prev/pop_prev were NA when the
  .h2_liab.log was read first); local modules aligned with nf-core specs ($args, full version reporting,
  nf-core-format meta.yml, dependency tags). `nf-core modules lint`: 0 failures.
- Full suite passes locally in Docker (nf-test 0.9.5, Nextflow 26.10.5): 23 module, 6 function,
  3 pipeline tests; pytest 48 passed / 6 skipped. CI still pins NFT_VER 0.9.4.
- Genome builds: all methods match variants by rsID and take positions from their references, so
  sumstats need rsIDs, not a particular build. Files without rsIDs get them from an rsID map
  (samplesheet `build`, `--rsid_map_grch37/38`; `.bim` or CHR/BP/SNP/A1/A2 table). No liftover.
- `tests/data/` is committed (`.gitignore` ignores only top-level `/data/`); it was missing on
  GitHub before 2026-11-09.
- Remaining nf-core gaps (warnings only): no meta.yml/tests for subworkflows method_ldsc/method_lava;
  nf-core utils subworkflows have updates; template TODOs; manifest.diagram; ro-crate.
- Pipeline not yet run in the published container on CI (image only local). No real genome-wide run yet.
- BigaGWAS's own LAVA scripts (specific_lava.R, 2traits.R, lava_pgc.R) are not public, so its settings are
  partly unknown.

## Possible next steps (only with the user's go-ahead)

1. Pelle/SLURM profile + first real genome-wide pair (e.g. MDD × another trait) — main next step for MapGEN.
   Expect ~1.5–3.5 CPU-hours per pair for LAVA.
2. CI workflow: build/push private gencorr-tools image to ghcr.io and run nf-test in it.
3. `lava/gencorr.lava_native.tsv` output in LAVA/BigaGWAS column layout.
4. Option to skip FRQ in LDSC munge; LAVA 0.1.0 container variant (exact BigaGWAS reproduction).
5. Another method (HDL or HDL-L) via docs/extending.md.
