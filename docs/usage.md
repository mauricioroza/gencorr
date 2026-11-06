# Usage

## 1. Samplesheet (`--input`)

CSV, one row per trait. Only `id` and `sumstats` are required.

| column | meaning |
|---|---|
| `id` | short trait ID: starts with a letter; letters, digits, `_ . -`. Becomes file names and the LAVA phenotype ID. |
| `sumstats` | path or URL to the summary statistics (plain or `.gz`; tab, comma or whitespace separated) |
| `n` | total sample size, used **only** if the file has no N column |
| `n_cases`, `n_controls` | binary traits. Used for N when missing, for LAVA's case/control model and for LDSC liability h2 |
| `pop_prev` | population prevalence (binary traits). Adds liability-scale h2 (LDSC) and LAVA `h2_latent` |
| `colmap` | manual column mapping, e.g. `Z=STAT,N=NMISS` (standard=source) |
| `label` | free-text display name (not used in computation) |

### Column detection

The harmonizer recognises the usual names case-insensitively (full list in
`bin/harmonize_sumstats.py`, `COLUMN_SYNONYMS`):

| standard | recognised source names (first match wins) |
|---|---|
| SNP | snp, rsid, rs_id, hm_rsid, variant_id, markername, snpid, snpid_ukb, … |
| A1 (effect) | a1, effect_allele, hm_effect_allele, ea, alt, allele1, tested_allele |
| A2 | a2, other_allele, hm_other_allele, nea, ref, allele2 |
| Z | z, zscore, z_score, z-score, z_stat |
| BETA | beta, hm_beta, b, effect, logor, log_or, log_odds |
| OR | or, odds_ratio, hm_odds_ratio (converted to BETA = ln OR) |
| SE, P, N, EAF, INFO | the usual spellings |

`STAT` is **never** mapped automatically. Depending on the software it is a Z,
t or chi-square statistic, so map it with `colmap` if you know it is a Z.

**A1 must be the effect allele.** GWAS Catalog harmonised files
(`hm_effect_allele`) are handled. Check files where "A1" means something else
(e.g. PLINK 1.9 `.assoc` A1 = minor allele, which *is* the tested allele).

### Signed Z-score

Per variant, the first available of (`--z_priority`, default `z,p,beta_se`):

1. an explicit Z column;
2. |Z| from the two-sided P, sign from BETA/ln OR;
3. BETA / SE.

Option 2 before 3 matches what LDSC `munge_sumstats.py` and LAVA do with raw
files (and what BigaGWAS fed them). That keeps results comparable.

## 2. Which pairs (`--pair_mode`)

| mode | pairs | BigaGWAS equivalent |
|---|---|---|
| `all` (default) | every pair of samplesheet traits | n/a |
| `target --targets mdd,scz` | each target vs every other trait | "massive" / one-vs-group jobs |
| `pairs --pairs pairs.csv` | exactly the listed pairs (`trait1,trait2` header) | "2 traits" jobs |

Every method analyses the same pairs.

## 3. Reference data

### LDSC

| file | where |
|---|---|
| European LD scores `eur_w_ld_chr/` | Zenodo record 8182036 (`eur_w_ld_chr.tar.gz`, copy of the Alkes group archive) |
| HapMap3 SNP list `w_hm3.snplist` | Zenodo record 7773502 (`w_hm3.snplist.gz`) |

```bash
--ldsc_ld /refs/eur_w_ld_chr          # files 1.l2.ldscore.gz ... 22.l2.ldscore.gz
--ldsc_merge_alleles /refs/w_hm3.snplist
# use other ancestries' LD scores for non-European GWAS (match the GWAS ancestry!)
```

`--ldsc_merge_alleles` is strongly recommended. Without it, LDSC fails a whole
pair as soon as one rsID has different alleles in the two traits.

### LAVA

| reference | notes |
|---|---|
| **UK Biobank, LAVA LD format v1.1** (recommended for European data) | links in the LAVA repository's `REFERENCE.md`; ~15 GB unzipped, per chromosome; requires LAVA ≥ 0.1.5 |
| 1000 Genomes PLINK files | smaller and quicker, but the LAVA authors report possible type I error inflation and biased local h2 with 1000G for European data |

```bash
--lava_ref /refs/lava-ukb/lava-ukb-v1.1    # prefix WITHOUT the _chrN part
--lava_ref /refs/g1000/g1000_eur           # PLINK prefix (.bed/.bim/.fam)
```

All files starting with the prefix are staged, so keep the reference in its own
directory.

Loci: the default `--lava_loci` is LAVA's 2,495 GRCh37 blocks (fetched from the
pinned LAVA commit). **Your summary statistics must be on the same build as the
reference and loci (GRCh37 for the defaults).** The pipeline does not lift over.

## 4. Important parameters

| parameter | default | |
|---|---|---|
| `--methods` | `ldsc,lava` | |
| `--lava_overlap` | `ldsc` | sample-overlap correction from LDSC intercepts (LAVA vignette recipe); `none`; or a matrix file |
| `--lava_snp_set` | `pairwise` | `pairwise`: a pair uses SNPs shared by its two traits. `joint`: all traits loaded together, SNPs shared by **all** traits (closer to a BigaGWAS group job) |
| `--lava_univ_thresh` | `0.05` | bivariate test only if both local h2 p-values are below this. For genome-wide screens the LAVA authors suggest a stricter (e.g. Bonferroni) filter |
| `--lava_loci_per_chunk` | `250` | parallelism only; results do not depend on it |
| `--lava_seed` | `20201231` | base seed for LAVA's Monte Carlo CIs/p-values |
| `--ldsc_rg_batch_size` | `50` | traits per `ldsc --rg` call; does not change estimates |
| `--drop_ambiguous` | `false` | remove A/T, C/G SNPs for all methods (LDSC munge always removes them) |
| `--min_info` | `0` (off) | remove variants with INFO ≤ value, for all methods (see 4a) |
| `--min_maf` | `0` (off) | remove variants with MAF ≤ value, for all methods (see 4a) |
| `--lava_snp_list` | none (all SNPs) | run LAVA only on the SNP IDs in this list, e.g. `w_hm3.snplist` (see 4b) |
| `--save_harmonized` / `--save_intermediates` | `false` | publish intermediate files |

`nextflow run mauricioroza/gencorr --help --show_hidden` lists everything.

## 4a. Variant QC: `--min_info`, `--min_maf`

Applied once, during harmonization, so **every method** sees the same
filtered variants. The rule is LDSC munge's: a variant is removed when its value
is **≤** the threshold. Missing values are kept (genotyped SNPs often have no
INFO).

* INFO comes from the INFO column (`info`, `info_ukb`, `impinfo`, `r2`, …).
* MAF comes from a `MAF` column if present, otherwise `min(EAF, 1 − EAF)`.
* If a file has no such column the threshold cannot be applied: nothing is
  removed and `harmonized/<id>.qc.json` contains `warning_min_info` /
  `warning_min_maf` (also printed in the task log).

**Recommended for real analyses: `--min_maf 0.01 --min_info 0.9`** (LDSC's
standard values). Without them, LDSC still filters MAF itself, but only when the
file has an allele-frequency column munge recognises; a file with only a
`MAF` column (like the vignette bmi file) would otherwise pass rare variants
into LDSC. LAVA applies no MAF/INFO filter of its own.

## 4b. SNP set for LAVA: `--lava_snp_list`

By default LAVA uses every harmonized SNP that is also in the LD reference.
`--lava_snp_list w_hm3.snplist` restricts LAVA to the IDs in a list
(first column; header optional; plain or gzipped). This is what BigaGWAS does
(HapMap3-based list), and it is needed to reproduce BigaGWAS' LAVA results.

| | all SNPs (default) | HapMap3 only (`--lava_snp_list w_hm3.snplist`) |
|---|---|---|
| information per locus | all SNPs; LAVA compresses LD to PCs itself | ~15–20% of SNPs; some loci drop to a handful of SNPs |
| robustness to imputation quality | depends on QC → use `--min_info`/`--min_maf` | common, well-imputed SNPs only |
| same SNPs across GWAS | pairwise intersection | fixed list |
| comparable with BigaGWAS | no | yes (shown on the test data: rho equal to ~1e-5) |
| ancestry | whatever the reference covers | HapMap3 is less complete outside European ancestry |

LDSC is not affected (it has its own list: `--ldsc_merge_alleles`). The SNP set
used is recorded in every LAVA row (`extra.snp_set`: `all` or the file name).

## 5. Running on UPPMAX (Pelle)

```bash
module load Nextflow           # or install Nextflow >= 25.10.4 yourself
export NXF_HOME=/proj/<project>/nobackup/$USER/.nextflow
export NXF_APPTAINER_CACHEDIR=/proj/<project>/nobackup/$USER/apptainer_cache

nextflow run mauricioroza/gencorr -profile apptainer -c pelle.config \
    --input samplesheet.csv --outdir results ...
```

with a small `pelle.config` (adjust account/partition to your allocation):

```groovy
process {
    executor       = 'slurm'
    clusterOptions = '-A <naiss_project>'
    resourceLimits = [ cpus: 16, memory: 120.GB, time: 240.h ]
}
executor { queueSize = 200; submitRateLimit = '10/1min' }
```

nf-core/configs also has an `uppmax` institutional profile (`-profile uppmax`),
loaded automatically from `--custom_config_base`. Check that it covers Pelle
before relying on it.

Memory: harmonization loads one GWAS file at a time in pandas (`process_low`,
12 GB). For very large files (> 20 M rows) raise it in a config:
`process { withName: 'HARMONIZE_SUMSTATS' { memory = 32.GB } }`.

`-resume` works as usual. Changing a LAVA parameter re-runs only LAVA.
