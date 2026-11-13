# LDSC on Python 3 (`--ldsc_impl py3`)

By default gencorr runs LDSC 1.0.1 (Bulik-Sullivan et al.), the last release of
the original code, which needs Python 2.7 and pandas 0.20 (bioconda package,
BioContainers image). `--ldsc_impl py3` runs the same steps (munge, h2, rg) with
[CBIIT/ldsc](https://github.com/CBIIT/ldsc), a port of that code to Python 3:

| | `original` (default) | `py3` |
|---|---|---|
| code | bulik/ldsc 1.0.1 | CBIIT/ldsc, branch `ldsc313`, commit `128f996` (2026-11-03) |
| logs say | `Version 1.0.1` | `Version 3.0.1` |
| Python, numpy, pandas, scipy | 2.7, 1.x, 0.20.3, 0.x | 3.13.7, 2.4.3, 2.3.3, 1.15.0 |
| image | `quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2` | `ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996` (`containers/gencorr-ldsc`, GPL-3.0) |
| conda | yes | no (container only) |

The port changes Python 2 idioms (`print`, `xrange`, `map` returning lists,
relative imports) and pandas/numpy API calls. The regression, jackknife and
munge code are the original ones. `LDSC_SUMMARISE` (the log parser) runs in
`gencorr-tools` with either implementation.

```bash
nextflow run . -profile docker --ldsc_impl py3 ...
```

## Tests inside the image

The image is built only if the port's own unit tests pass (`test/`, 167 tests,
run with pytest + [pynose](https://pypi.org/project/pynose/) because they import
`nose.tools`, and the original nose does not import on Python ≥ 3.12):

```bash
docker build -t ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996 containers/gencorr-ldsc
docker run --rm ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996 run_upstream_tests.sh
# 164 passed, 3 deselected
```

`run_upstream_tests.sh` changes only the test harness, never LDSC's code:
a nose-style `setUpClass` on a plain class renamed to pytest's `setup_class`;
the simulated LD-score fixtures gzipped (the port reads only `.l2.ldscore.gz`);
one test that passed a 1×1 `np.matrix` where LDSC passes a float. The 3
deselected tests check behaviour the port changed (see below).

gencorr's own Python tests also pass inside this image, including the LDSC
oracle tests (`tests/python/test_ldsc_oracle.py`, which then run the port):

```bash
docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp -v $PWD:/w -w /w \
    ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996 python -m pytest -q tests/python
# 93 passed, 6 skipped (LAVA and nextflow tests: no R/nextflow in this image)
```

## Comparison with LDSC 1.0.1 and BigaGWAS

The same inputs were run with both implementations (`--save_intermediates true`)
and compared output by output with `tests/bigagwas/compare_ldsc_impl.py`
(results table incl. every `extra` field, h2 table, sample-overlap matrix,
munged files row by row, munge SNP counts). References as in
[bigagwas_comparison.md](bigagwas_comparison.md): `eur_w_ld_chr`, `w_hm3.snplist`.

| data | pairs | result |
|---|---|---|
| BigaGWAS job 971 files (T2D × CHD), as for job 972 | 1 | **identical**: munged Z and N for all 1,217,311 rows of both traits, rg, se, z, p, h2, intercepts, gencov intercept, 2,937 SNPs, overlap matrix |
| same, vs BigaGWAS job 972's own LDSC log | 1 | **identical**: rg 0.3476, se 0.0796, z 4.3684, p 1.2518e-05, 2,937 SNPs, h2 intercepts 1.7697 / 0.474, gencov intercept 0.0643 |
| BigaGWAS upload files (bmi, neuro, depression, rheuma), incl. job 966 (bmi × depression, rg 0.1154) | 6 | munged files, rg, se, z, h2 table and overlap matrix identical; p identical in 5/6 pairs, neuro × rheuma 0.0375 vs 0.037473 (printing, see below) |
| `-profile test` (small test references) | 6 | rg, se, z, h2, overlap and every LAVA estimate identical; p as above |

Run time and memory on these (small) inputs: h2 and rg the same (4–7 s), munge
~25% faster with the port (7.5 s vs 10 s), peak memory ~10–30% higher
(~535 MB vs 400–500 MB). Genome-wide inputs have not been timed yet.

Reproduce (WSL, references in `~/refs`):

```bash
for impl in original py3; do
  nextflow run . -profile docker -c local.config --input samplesheet_971.csv \
    --outdir biga_ldsc_impl/971_$impl --methods ldsc \
    --ldsc_ld ~/refs/eur_w_ld_chr --ldsc_merge_alleles ~/refs/w_hm3.snplist \
    --ldsc_impl $impl --save_intermediates true
done
tests/bigagwas/compare_ldsc_impl.py --a biga_ldsc_impl/971_original --b biga_ldsc_impl/971_py3
```

## Known differences of the port

Found by running both versions and by reading the diff against LDSC 1.0.1
(bulik/ldsc `0a7ebd9`). Only the first affects gencorr's numbers, and
`parse_ldsc.py` handles it.

1. **p in the rg summary table.** The table is printed with pandas'
   `DataFrame.to_string`. pandas 0.20 printed the `p` column in scientific
   notation (`1.2518e-05`, `5.1006e-02`); current pandas prints 4 fixed
   decimals, so p < 5e-5 appears as `0.0000`. `parse_ldsc.py` then takes p from
   the pair's own `P:` line, which both versions print with 4 significant
   digits (`P: 1.2518e-05`). Result: small p are exact; larger p have 4 instead
   of 5 significant digits (0.0375 vs 0.037473).
2. **LD score files must be gzipped.** `parse.ldscore` looks only for
   `*.l2.ldscore.gz`; plain or `.bz2` LD scores (accepted by 1.0.1) are not
   found. The published `eur_w_ld_chr`/`eas_ldscores` files are `.gz`.
3. **Missing chromosomes are skipped.** With `--ref-ld-chr`/`--w-ld-chr` the
   port uses the chromosomes whose files exist; 1.0.1 stops when one of 1–22 is
   missing. Check `ls <ldsc_ld>` if a reference set might be incomplete.
4. **rg when a jackknife denominator is negative.** 1.0.1 raises a
   floating-point error (gencorr status `ldsc_error`); the port computes it
   under `np.errstate(invalid='ignore')`, so rg is reported with a missing SE
   and p (gencorr status `ok`, `rg_se` NA). This happens only when a trait's
   h2 is close to zero.
5. **Unreadable `.gz` input files.** The port prints the pandas error and
   continues, so LDSC fails later with a less clear message. (Deselected test
   `test_read_sumstats`.)
6. Cosmetic: floats in munge logs are printed with Python 3's 17 digits
   (`N < 71533.33333333333`), and pandas prints a FutureWarning from
   `allele_merge`; neither changes the output files.
