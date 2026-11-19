# Multivariate GWAMA (`--gwama_groups`)

gencorr can combine a group of genetically correlated traits into **one new
trait** with N-weighted multivariate GWAMA (N-GWAMA; Baselmans et al. 2019,
[doi:10.1038/s41588-018-0320-8](https://doi.org/10.1038/s41588-018-0320-8)).
The combined trait has the same harmonized format as every other trait, so
LDSC, LAVA, SumHer and Popcorn analyse it like a samplesheet trait. For
example, combine life satisfaction, positive affect, neuroticism and
depressive symptoms into a well-being trait and estimate its rg with other
traits.

GWAMA does not estimate rg itself. It makes new traits, before the methods
run.

## Usage

```csv title="groups.csv"
group,trait,sign
wellbeing,ls,1
wellbeing,pa,1
wellbeing,neu,-1
wellbeing,dep,-1
```

```bash
nextflow run mauricioroza/gencorr -profile docker \
    --input samplesheet.csv --outdir results \
    --gwama_groups groups.csv \
    --methods ldsc,lava --ldsc_ld eur_w_ld_chr/ --ldsc_merge_alleles w_hm3.snplist \
    --lava_ref g1000_eur
```

- `group`: name of the combined trait. It follows the samplesheet ID rules and
  must not be a samplesheet trait ID. A file can have several groups, and a
  trait can be in several groups.
- `trait`: a samplesheet trait ID. A group needs at least 2 traits, all of one
  ancestry.
- `sign` (optional, default `1`): `-1` reverses the trait's effects before
  combining. The traits of a group should correlate positively: Baselmans et
  al. reversed neuroticism and depressive symptoms to combine them with life
  satisfaction and positive affect.

The combined traits join the pairs like samplesheet traits: with
`--pair_mode all` every trait is paired with every combined trait, and
`--targets` and `--pairs` may name combined traits. A combined trait paired
with one of its own traits gives a high rg by construction.

`--gwama_groups` needs LDSC LD scores (`--ldsc_ld`, or `--ldsc_ld_eas` for an
EAS group) even if `ldsc` is not in `--methods`, and a container profile
(LDSC has no conda package).

| parameter        | default     |                                                                                |
| ---------------- | ----------- | ------------------------------------------------------------------------------ |
| `--gwama_groups` | none        | CSV with `group,trait[,sign]`                                                  |
| `--gwama_n`      | `effective` | N column of a combined trait: `effective` or `observed` (sum of the traits' N) |
| `--gwama_min_rg` | `0.5`       | warn when two traits of a group have an LDSC rg (after signs) below this       |

## What the pipeline does

1. **LDSC on the group's traits.** Every pair inside each group goes through
   LDSC (munge, h2, rg; one LDSC run for all groups). This gives each trait's
   observed-scale SNP heritability h2 and univariate intercept, and each pair's
   cross-trait intercept `gcov_int` and rg. Outputs: `gwama/ldsc/`.
2. **N-GWAMA per group** (`GWAMA_COMBINE`, `bin/gwama.py`). Per SNP, over the
   traits that have it, with weights w_i = N_i · h2_i:

   ```
   Z    = Σ_i √w_i Z_i  /  √( Σ_ij √w_i √w_j C_ij )
   N    = √N' C⁻¹ √N                          (--gwama_n effective)
        = Σ_i N_i                              (--gwama_n observed)
   EAF  = N-weighted mean of the traits' EAF
   BETA = Z · SE,   SE = 1 / √(N_obs · 2 EAF (1 − EAF))
   P    = two-sided, from Z
   ```

   C is the cross-trait intercept (CTI) matrix: the univariate LDSC intercepts
   on the diagonal and `gcov_int` off the diagonal (averaged when LDSC ran a
   pair in both orientations). It corrects for sample overlap between the
   GWAS and for inflation from population stratification. A trait with sign
   −1 has its Z-scores, its `gcov_int` and its rg with the others reversed.

3. **The combined trait joins the harmonized traits** and every method runs
   on it.

SNPs are combined over the union of the traits: a SNP that is missing in some
traits is combined from the others (its N and effective N are lower). Alleles
are aligned to the first trait in the group that has the SNP; a SNP whose
alleles do not match that trait's, not even swapped, is left out for that
trait only.

## Outputs

```
results/gwama/
├── <group>.harmonized.tsv.gz   the combined trait (harmonized format; INFO is NA)
├── <group>.gwama.json          inputs and QC (below)
└── ldsc/                       LDSC run on the groups' pairs
    ├── gencorr.ldsc.tsv        rg and gcov_int of every pair inside a group
    ├── gencorr.ldsc_h2.tsv     h2 and intercept of every trait in a group
    ├── gencorr.sample_overlap.txt
    ├── logs/
    └── munged/
```

`<group>.gwama.json` records, per trait, the sign, the h2 used and the
alignment counts (`aligned`, `flipped`, `allele_mismatch`); the CTI and rg
matrices used (after signs); `rows_out`; `snps_in_all_traits` and
`snps_by_n_traits`; `genome_wide_significant` (P ≤ 5e-8); and `warnings`.

## When it helps, and warnings

Baselmans et al.'s simulations: N-GWAMA gains power over single-trait GWAMA
when the traits' rg is about 0.5 or higher, and the gain disappears at
rg ≤ 0.4. When a SNP affects fewer than half of the traits, single-trait
analysis can find more. The pipeline warns (task log and `warnings` in the
JSON) when a pair in a group has an rg below `--gwama_min_rg`, a negative rg
after signs (give one of the traits sign `-1`), or no rg estimate.

The pipeline stops when a trait's h2 is not positive (N-GWAMA weights traits by
N · h2), when the LDSC intercept is missing (LDSC run with a constrained
intercept), when a pair's `gcov_int` is missing (e.g. the LDSC rg failed) or
when the CTI matrix is not positive definite.

A combined trait has no case/control counts or prevalence: methods treat it as
a continuous trait (LDSC reports observed-scale h2 only).

## Validation

`tests/python/test_gwama.py` compares `bin/gwama.py` with the authors' R
function (`N_weighted_GWAMA.function.1_2_6.R`,
[baselmans/multivariate_GWAMA](https://github.com/baselmans/multivariate_GWAMA)
at commit `d55b617`) on their example data (4 traits, 100,000 SNPs, their CTI
matrix and h2 values). Z, P, effective N, BETA, SE and EAF agree to a relative
difference below 1e-6 for all SNPs. The R function and the example data are
downloaded at test time (the repository has no licence, so nothing from it is
copied into gencorr); the test runs where Rscript has data.table, e.g. in the
`gencorr-tools` image, and is skipped elsewhere.

Known differences from the R function:

- **Missing EAF.** The R function drops every input row with a missing value,
  so a SNP without EAF in one trait is combined from the other traits only.
  gencorr's harmonized format allows a missing EAF, so `gwama.py` keeps the
  SNP's Z and N and averages EAF over the traits that have it (BETA and SE are
  NA when no trait has an EAF). The oracle test removes those rows (45 in three
  of the example files) from both inputs.
- **Not included:** the R function's options `method='intersection'`,
  `covars` and `type='demontis'` (mixing case-control and continuous traits),
  and MA-GWAMA, the paper's model-averaging method (one `metafor` model per
  SNP and candidate model; too slow genome-wide, and the authors say LDSC
  results on its output should be read with caution).
- **CTI from the pipeline's LDSC.** Their example uses a CTI matrix and h2
  values computed beforehand; gencorr computes both with LDSC on the same
  harmonized files the methods use.
