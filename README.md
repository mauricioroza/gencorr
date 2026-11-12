# mauricioroza/gencorr

Global and local **genetic correlation** from GWAS summary statistics: a
reproducible Nextflow pipeline (nf-core template) that runs LDSC and LAVA on any
set of traits and returns one tidy results table.

```
   __ _  ___ _ __   ___ ___  _ __ _ __
  / _` |/ _ \ '_ \ / __/ _ \| '__| '__|
 | (_| |  __/ | | | (_| (_) | |  | |
  \__, |\___|_| |_|\___\___/|_|  |_|
  |___/   global + local genetic correlation
```

| | method | tool | what you get |
|---|---|---|---|
| **global** | genome-wide rg | [LDSC](https://github.com/bulik/ldsc) 1.0.1 (same version as BigaGWAS) | rg, SE, p per trait pair; observed- and liability-scale h2 per trait |
| **global** | genome-wide rg, LDAK heritability model | [SumHer](https://dougspeed.com/genetic-correlations/) (LDAK 6.1) | rg, SE, p per trait pair; h2 and coheritability |
| **cross-ancestry** | rg between a EUR and an EAS GWAS | [Popcorn](https://github.com/brielin/popcorn) (commit `20a29e8`) | rg (genetic impact or effect), SE, p; h2 in each population |
| **local** | rg per genomic locus | [LAVA](https://github.com/josefin-werme/LAVA) 0.1.5 (commit `e729a24`) | local rho, 95% CI, p per pair and locus; local h2 per trait |

It was built to replace [BigaGWAS](https://bigagwas.org) for routine analyses:

* **Same numbers as BigaGWAS.** With the same data and settings, LDSC rg is
  identical and LAVA rho agrees to 1e-5–1e-7
  ([comparison](docs/bigagwas_comparison.md#6-results-on-the-test-jobs)).
* **Faster.** No rsID mapping step, LAVA loci run in parallel chunks and empty
  loci are skipped. On the example data one LDSC job took ~45 s instead of
  ~20 min on BigaGWAS.
* **Reproducible.** Every tool version is pinned, LAVA's Monte Carlo steps are
  seeded per (locus, pair), and outputs are byte-identical on re-runs.
* **Nothing disappears.** Pairs and loci that could not be estimated stay in
  the table with a `status` explaining why.
* **Extensible.** Every method reads the same harmonized input and writes the
  same table, so adding HDL, SumHer, Popcorn or HDL-L is one subworkflow
  ([docs/extending.md](docs/extending.md)).

## How it works

![gencorr metro map: GWAS summary statistics (with an optional rsID map) are harmonized, then same-ancestry pairs run through LDSC and SumHer (global rg) and LAVA (local rg, corrected with LDSC's sample-overlap matrix), EUR × EAS pairs through Popcorn (cross-ancestry rg), and everything is merged into gencorr_results.tsv](docs/images/gencorr_metro_map.svg)

1. **Harmonize** each GWAS once: column names are detected automatically
   (SNP/rsid/MarkerName, A1/effect_allele, BETA/OR/LogOR, SE, P, Z, N, …),
   OR becomes ln OR, a signed Z is derived, and duplicates, non-SNVs and
   unusable rows are removed. A QC report per trait records every decision.
2. **LDSC** munges, estimates h2 per trait and rg for every pair, and builds
   the sample-overlap matrix from its intercepts.
3. **LAVA** tests every pair in each of LAVA's 2,495 GRCh37 loci, corrected for
   sample overlap by default.

   *Sample overlap* means the original GWAS shared participants (e.g. UK Biobank
   in both). Their estimation errors are then correlated, which inflates rg.
   LDSC's cross-trait intercepts measure this from the summary statistics alone
   (no individual-level data needed); LDSC corrects itself, and gencorr passes
   the intercepts to LAVA as its overlap matrix (`--lava_overlap ldsc`; `none`
   to switch off).
4. **Merge** checks every table against one schema and writes
   `gencorr_results.tsv`.

## Quick start

**1. Requirements:** [Nextflow](https://www.nextflow.io/) ≥ 25.10.4 and Docker,
Apptainer or Singularity.

**2. Build the helper image once.** LAVA's licence (*All rights reserved*) does
not allow a public image, so `ghcr.io/mauricioroza/gencorr-tools:0.1.0` is
private. Build it from [containers/gencorr-tools](containers/gencorr-tools/Dockerfile),
or log in to ghcr.io with access to it:

```bash
docker build -t ghcr.io/mauricioroza/gencorr-tools:0.1.0 containers/gencorr-tools
```

LDSC comes from the public BioContainers image and needs no setup. The
SumHer and Popcorn images (`ghcr.io/mauricioroza/gencorr-ldak:6.1`,
`ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8`, both GPL-3.0) are
published on ghcr.io and pulled automatically; their recipes are in
[containers/](containers/).

**3. Check your setup** on the bundled test data (4 traits, 6 pairs, a few minutes):

```bash
nextflow run mauricioroza/gencorr -profile test,docker --outdir test_out
```

**4. Real run** (use `-profile apptainer` on a cluster):

```bash
nextflow run mauricioroza/gencorr -profile docker \
    --input samplesheet.csv --outdir results \
    --ldsc_ld /refs/eur_w_ld_chr --ldsc_merge_alleles /refs/w_hm3.snplist \
    --lava_ref /refs/lava-ukb/lava-ukb-v1.1 \
    --min_maf 0.01 --min_info 0.9
```

Where to download the LD scores, SNP list and LAVA reference:
[docs/usage.md § Reference data](docs/usage.md#3-reference-data).

### Samplesheet

One row per trait; only `id` and `sumstats` are required:

```csv
id,sumstats,n,n_cases,n_controls,pop_prev,colmap,label
bmi,/data/bmi.txt.gz,,,,,,Body mass index
mdd,/data/mdd.tsv.gz,,170756,329443,0.15,,Major depression
scz,/data/scz.txt,,53386,77258,0.01,"Z=STAT",Schizophrenia
t2d,GCST000028,,,,,,Type 2 diabetes (GWAS Catalog)
```

* `sumstats` is a path, a URL, or a **GWAS Catalog accession** (`GCST…`): the
  study's harmonised file is downloaded and md5-checked, and its sample size
  is used when `n` is empty and the study has one sample group, of the trait's ancestry
  ([details](docs/usage.md#gwas-catalog-studies)).
* `n` is used only when the file has no N column.
* `n_cases`/`n_controls` mark a binary trait (LAVA case/control model, N fallback);
  with `pop_prev` you also get liability-scale h2.
* `colmap` overrides column detection, e.g. `Z=STAT`. `STAT` is never mapped
  automatically, because it is not always a Z-score.

### Common options

| option | default | |
|---|---|---|
| `--methods` | `ldsc,lava` | run one or both. With `--methods lava` alone, also set `--lava_overlap none` (the default overlap comes from LDSC) |
| `--pair_mode` | `all` | `all` pairs; `target --targets mdd,scz` = each target vs all others (BigaGWAS "one vs group"); `pairs --pairs pairs.csv` = listed pairs only |
| `--min_maf`, `--min_info` | off | variant QC for all methods; `0.01` / `0.9` recommended |
| `--lava_snp_list` | all SNPs | e.g. `w_hm3.snplist`, to reproduce BigaGWAS' LAVA |
| `--lava_univ_thresh` | `0.05` | local h2 p-value both traits must pass before a pair is tested |

All parameters: [docs/usage.md](docs/usage.md) or
`nextflow run mauricioroza/gencorr --help --show_hidden`.

## Outputs

```
results/
├── gencorr_results.tsv     all methods, one row per (method, trait pair, locus)
├── harmonized/             QC report per trait
├── ldsc/                   per-trait h2, sample-overlap matrix, LDSC logs
├── lava/                   LAVA rows, local h2 per locus and trait
└── pipeline_info/          exact software versions, parameters, run reports
```

`gencorr_results.tsv`, main columns (all columns: [docs/output.md](docs/output.md)):

| method | scope | trait1 | trait2 | locus | rg | rg_se | rg_ci_low | rg_ci_high | p | status | extra |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ldsc | global | bmi | mdd | NA | 0.11 | 0.07 | … | … | 0.11 | ok | {"gcov_int": …} |
| lava | local | bmi | mdd | 464 | 0.51 | NA | 0.27 | 0.81 | 1.1e-4 | ok | {"r2": …} |
| lava | local | bmi | mdd | 100 | NA | NA | NA | NA | NA | univ_filtered | … |

`status` is `ok` when there is an estimate; otherwise it says why not
(`univ_filtered`, `no_snps_in_locus`, `ldsc_error`, `estimate_out_of_bounds`, …).
Read the file in R with `read.delim("gencorr_results.tsv", quote = "")` and the
JSON `extra` column with `jsonlite::fromJSON`.

## Good to know

* **Genome build:** variants are matched by rsID, so summary statistics on
  GRCh37 or GRCh38 work without liftover if they have rsIDs. Files without
  rsIDs (`1:12345:A:G`-style IDs) get them from an rsID map for their build:
  set `build` in the samplesheet and `--rsid_map_grch37` / `--rsid_map_grch38`
  ([details](docs/usage.md#genome-build-and-rsids)).
* **Ancestry:** set `ancestry` (`EUR` default, `EAS`) per trait in the
  samplesheet. Each same-ancestry pair runs with that ancestry's references
  (`--ldsc_ld`/`--lava_ref` for EUR, `--ldsc_ld_eas`/`--lava_ref_eas` for EAS);
  EUR × EAS pairs are estimated by Popcorn (`--methods popcorn`) and reported
  by the other methods with status `cross_ancestry`
  ([details](docs/usage.md#ancestry)). For European data the LAVA authors
  recommend their UK Biobank reference over 1000 Genomes.
* **conda:** LAVA has no conda package, so `-profile conda` works only with
  `--methods ldsc`.
* **Clusters:** see [docs/usage.md § UPPMAX](docs/usage.md#5-running-on-uppmax-pelle)
  for a SLURM config. `-resume` works as usual; changing a LAVA option re-runs
  only LAVA.
* **Status:** v0.1.0dev, not yet released. Validated on the LAVA vignette and
  BigaGWAS example data; no genome-wide production run yet.

## Testing

54 Python tests, 23 nf-test module tests, 6 function tests and 3 pipeline
tests, including *oracle* tests that run LDSC and LAVA directly and check that
the pipeline gives the same numbers. `nf-core pipelines lint` and
`nf-core modules lint` report no failures. How to run them:
[docs/validation.md](docs/validation.md).

## Documentation

| | |
|---|---|
| [usage.md](docs/usage.md) | samplesheet, all parameters, reference data, running on UPPMAX |
| [output.md](docs/output.md) | every output file and column, status values |
| [design.md](docs/design.md) | design decisions, and where results differ from BigaGWAS on purpose |
| [validation.md](docs/validation.md) | test layers and how to run them |
| [bigagwas_comparison.md](docs/bigagwas_comparison.md) | reproducing BigaGWAS jobs and the comparison results |
| [extending.md](docs/extending.md) | adding a new method, step by step |
| [CONTRIBUTING.md](docs/CONTRIBUTING.md) | conventions and checks before a pull request |

## Credits and citation

Written by MR (Genis). The method choices follow BigaGWAS, whose LDSC and LAVA
usage this pipeline reproduces and extends.

If you use gencorr, please cite the tools it runs: LDSC
(Bulik-Sullivan et al. 2015, [10.1038/ng.3211](https://doi.org/10.1038/ng.3211);
[10.1038/ng.3406](https://doi.org/10.1038/ng.3406)), LAVA (Werme et al. 2022,
[10.1038/s41588-022-01017-y](https://doi.org/10.1038/s41588-022-01017-y)),
Nextflow and nf-core. Full list: [CITATIONS.md](CITATIONS.md).
