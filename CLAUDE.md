# gencorr — notes for Claude

Nextflow (nf-core template) pipeline for global and local genetic correlation (LDSC + LAVA)
for the Genis group. Owner: GitHub `mauricioroza`,
repo `mauricioroza/gencorr` (private), working branch `dev` (also `master`, `TEMPLATE`).

Goals: reproducible, testable, well documented; results checkable against the original tools;
new methods easy to add. A Rust rewrite was considered and rejected (stats run inside LDSC/LAVA;
the slow part of such platforms is rsID mapping, which this pipeline avoids by design).

Docs/code framing: BigaGWAS (BIGA) is mentioned ONLY in the acknowledgement (README Credits,
CITATIONS "Related work") and where its results were used to compare/test gencorr
(docs/bigagwas_comparison.md incl. its "Known differences in defaults" table, docs/validation.md
comparison sections, tests/bigagwas/, compare_results.py `bigagwas_ldsc` format). Don't describe gencorr
as replacing or following it anywhere else (docs, help texts, code comments). The notes below
mention it only as a comparison reference.

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
- Boolean params: with the default (v2) parser `--x false` on the CLI arrives as the String "false"
  (truthy); under NXF_SYNTAX_PARSER=v1 (nf-test) it is a real Boolean. Always test
  `params.x.toString() == 'true'`, never `params.x ? ...`. Guarded by tests/python/test_cli_options.py
  (`GENCORR_NEXTFLOW_TESTS=1`). Template's monochrome_logs still has this (cosmetic).
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
  (gzip mtime 0, `%.8g`) — preserve this. Memory: load_input reads only mapped columns (numerics
  parsed while reading, CHR/alleles shared strings), filters use one keep mask; 9.5M rows ~4.4 GB.
- LDSC 1.0.1 (bioconda, same as BigaGWAS): munge with merge-alleles (`--signed-sumstats Z,0`,
  `--chunksize 500000`: the container's pandas 0.20.3 `isin` on >1M rows uses np.in1d, O(n*m) for
  rsIDs, hours on genome-wide files; keep chunks <= 1M), h2
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
- Python unit + oracle: 99 tests (`tests/python`; the oracle tests skip unless LDSC/LAVA are on PATH).
  Run with the venv active: `source ~/venvs/gencorr/bin/activate && pytest tests/python` (scripts are
  called as `python3`, so the venv must be first on PATH). nf-test: module (34), function (9),
  pipeline (6).
- nf-core tools 4.2.1 is in `~/venvs/gencorr`: `nf-core pipelines lint`, `nf-core modules lint --all --local`
  (run from the repo root). Module warnings about the private container, custom meta keys, underscores in
  names and `process_single`+`process_long` are expected for this private pipeline.
- Test data is downloaded from raw.githubusercontent.com on every pipeline test; GitHub can throttle it
  (seen 2026-11-07, ~150 KB/s), which makes Nextflow hang on staging. Check speed with curl first.
- nf-test snapshots must contain only machine-independent values (identical under pandas 2.2.3 and 3.0.5).
- SumHer vs BigaGWAS job 974 (2traits SumHer on job 971's T2D x CHD, LDAK 5.2, CHR:BP tagging file,
  --genomic-control YES): gencorr (LDAK 6.1, rsID tagging file) identical to all 6 printed decimals.
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
- Full suite passes locally in Docker (nf-test 0.9.5, Nextflow 26.10.5): 31 module, 6 function,
  5 pipeline tests; pytest 60 passed / 6 skipped. CI still pins NFT_VER 0.9.4.
- Popcorn (`--methods popcorn`, cross-ancestry pairs only; registry `cross_ancestry: true`): image
  `ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8` (GPL-3.0; pushed to ghcr.io 2026-11-09). Needs
  numpy<1.24 (np.bool) and six (pysnptools). Popcorn's p tests rg<1 -> kept in extra; `p` = rg!=0.
  Container-only (conda_ok false). Test data: Popcorn repo test files at the pinned commit + score file
  tests/data/popcorn. `.git/info/exclude` hides samplesheet_*.csv: add test samplesheets with `git add -f`.
- SumHer (LDAK 6.1, `--methods sumher`): image `ghcr.io/mauricioroza/gencorr-ldak:6.1`
  from containers/gencorr-ldak (GPL-3.0; pushed to ghcr.io 2026-11-09). LDAK 6 writes `SE` not `SD`,
  can't read gzipped tagging files (module gunzips), defaults to multiplicative-inflation model (we pass
  `--genomic-control NO` like BigaGWAS). LDAK's current pre-computed tagging files (downloaded 2026-11-09)
  use rsIDs; older ones (BigaGWAS) use GRCh37 CHR:BP -> `--sumher_snp_positions` only for those. `--calc-tagging` needs ~9 GB RAM by default: use `--bit-size 400` in WSL.
- Ancestry: samplesheet `ancestry` (EUR default, EAS); refs `--ldsc_ld[_eas]`, `--ldsc_w_ld[_eas]`,
  `--lava_ref[_eas]` (see ancestryParam/traitAncestry in utils subworkflow). Same-ancestry pairs only go
  to methods; cross-ancestry pairs become `cross_ancestry` rows (gencorr.nf notEstimatedRow).
- Genome builds: all methods match variants by rsID and take positions from their references, so
  sumstats need rsIDs, not a particular build. Files without rsIDs get them from an rsID map
  (samplesheet `build`, `--rsid_map_grch37/38`; `.bim` or CHR/BP/SNP/A1/A2 table). No liftover.
- GWAS Catalog accessions in samplesheet `sumstats` (GCST\d+ -> meta.gwas_catalog, resolveSumstats in utils):
  GWASCATALOG_FETCH (bin/fetch_gwascatalog.py, stdlib only incl. a small reader for the Catalog's
  -meta.yaml) -> HARMONIZE_SUMSTATS 4th input source_info JSON (provenance in qc.json `source`; N fallback
  only when the metadata lists ONE sample group of the samplesheet ancestry; real studies list the same
  sample twice (GCST90002409) or 'NR' ancestry (GCST004131), so sums are never used). Old hm_ files: hm_ alleles+effects
  used together. `--gwascatalog_cache` (storeDir), `--gwascatalog_base` (mirror/local folder; tests use
  synthetic tests/data/gwascatalog). BigaGWAS's OpenGWAS route is dead (old URLs 404, API needs JWT
  since 2024-05); its Neale/FinnGen/PGC sources are server-local copies.
- `tests/data/` is committed (`.gitignore` ignores only top-level `/data/`); it was missing on
  GitHub before 2026-11-09.
- Remaining nf-core gaps (warnings only): no meta.yml/tests for subworkflows method_ldsc/method_lava;
  nf-core utils subworkflows have updates; template TODOs; manifest.diagram; ro-crate.
- Images on ghcr.io (2026-11-09): gencorr-ldak:6.1 (47 MB compressed) and gencorr-popcorn:1.0.0-20a29e8
  (225 MB); both were created private (GitHub default) and may be made public by the user (GPL-3.0).
  gencorr-tools:0.1.0 is NOT on ghcr.io (only local) and must stay private when pushed (LAVA).
  Pushing needs `gh auth refresh -s write:packages` + `gh auth token | docker login ghcr.io -u mauricioroza
  --password-stdin`, run by the user in the Ubuntu terminal (not sudo, not PowerShell).
- Branch `ldsc-py3` (2026-11-13): `--ldsc_impl py3` = CBIIT/ldsc branch ldsc313 @128f996 (Python 3.13,
  numpy 2.4.3, pandas 2.3.3) in `containers/gencorr-ldsc`, image `ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996`
  (GPL-3.0; built locally, NOT pushed yet). Swapped in via `conf/ldsc_py3.config` (conditional includeConfig;
  strict config forbids `if` blocks). Upstream tests run at image build (`run_upstream_tests.sh`: 164 passed,
  3 deselected). Identical to 1.0.1 on BigaGWAS data (job 972 exact) except p printing: new pandas rounds the
  summary-table p to 4 decimals, so parse_ldsc.py reads the pair's `P:` line. Details: docs/ldsc_py3.md.
- Pipeline not yet run in the published containers on CI. No real genome-wide run yet.
- BigaGWAS's own LAVA scripts (specific_lava.R, 2traits.R, lava_pgc.R) are not public, so its settings are
  partly unknown.

## Possible next steps (only with the user's go-ahead)

1. First real genome-wide pair (e.g. MDD × another trait). Expect ~1.5–3.5 CPU-hours per pair for LAVA.
   Remote compute is still to be found (Pelle/UPPMAX was the user's old SU job: no longer available;
   docs/usage.md §5 is a generic SLURM example).
2. CI workflow: build/push private gencorr-tools image to ghcr.io and run nf-test in it.
3. `lava/gencorr.lava_native.tsv` output in LAVA/BigaGWAS column layout.
4. Option to skip FRQ in LDSC munge; LAVA 0.1.0 container variant (exact BigaGWAS reproduction).
5. Another method (HDL or HDL-L) via docs/extending.md.
