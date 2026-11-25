# Usage

## 1. Samplesheet (`--input`)

CSV, one row per trait. Only `id` and `sumstats` are required.

| column                  | meaning                                                                                                                                                                                                |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `id`                    | short trait ID: starts with a letter; letters, digits, `_ . -`; no `__` (it joins IDs in pair names) and no trailing `_`. Becomes file names and the LAVA phenotype ID.                                |
| `sumstats`              | path or URL to the summary statistics (plain or `.gz`; tab, comma or whitespace separated), **or a GWAS Catalog accession** such as `GCST90000016` (see [GWAS Catalog studies](#gwas-catalog-studies)) |
| `n`                     | total sample size, used **only** if the file has no N column                                                                                                                                           |
| `n_cases`, `n_controls` | binary traits. Used for N when missing, for LAVA's case/control model and for LDSC liability h2                                                                                                        |
| `pop_prev`              | population prevalence (binary traits). Adds liability-scale h2 (LDSC) and LAVA `h2_latent`                                                                                                             |
| `colmap`                | manual column mapping, e.g. `Z=STAT,N=NMISS` (standard=source)                                                                                                                                         |
| `ancestry`              | `EUR` (default), `EAS`, `AFR`, `SAS`, `AMR` or `MULTI`: ancestry of the GWAS. Selects the LD references and decides which pairs each method analyses (see [Ancestry](#ancestry))                       |
| `build`                 | `GRCh37` (default) or `GRCh38`: build of the file's positions. Only used to pick the rsID map for variants without an rsID (see [Genome build and rsIDs](#genome-build-and-rsids))                     |
| `label`                 | free-text display name (not used in computation)                                                                                                                                                       |

### Ancestry

LDSC and LAVA model the LD between nearby variants, and LD differs between
populations, so each needs an LD reference of the **same ancestry as the GWAS**.
They also assume **both traits of a pair come from one population**.

| pair                                  | what happens                                                                                                                         |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| EUR × EUR                             | LDSC with `--ldsc_ld`, LAVA with `--lava_ref`, SumHer with `--sumher_tagfile`                                                        |
| X × X (X = EAS, AFR, SAS, AMR)        | LDSC with `--ldsc_ld_<x>` (`--ldsc_w_ld_<x>`), LAVA with `--lava_ref_<x>`, SumHer with `--sumher_tagfile_<x>` (e.g. `--ldsc_ld_afr`) |
| different ancestries (e.g. EUR × EAS) | **Popcorn** (`--methods popcorn`, see below); LDSC, LAVA and SumHer report the pair with `status = cross_ancestry`                   |
| any pair with a `MULTI` trait         | not estimated: every method reports `status = multi_ancestry`                                                                        |

A reference is only required if some pair needs it: a run where every trait is
EAS needs only the `_eas` references. `--ldsc_merge_alleles` (HapMap3) and
`--lava_loci` are shared by all ancestries; the default LAVA loci were defined
from European LD (GRCh37) and are commonly reused for other ancestries, which
keeps local results comparable across ancestries (block boundaries are not
LD-independent in AFR). With `--lava_snp_set joint`, traits are loaded together
per ancestry.

**`MULTI`** is for multi-ancestry meta-analyses (e.g. the multi-ancestry PGC
MDD results). No LD reference matches a mixed sample, so LDSC, LAVA, SumHer and
Popcorn would give biased estimates; gencorr harmonizes the file (QC report)
but analyses none of its pairs. Use the ancestry-specific results instead.

**`AMR`** (admixed American, e.g. Hispanic/Latino cohorts): LD scores from an
external reference underestimate h2 by 10–60% in admixed samples (Luo et al.
2021). Prefer covariate-adjusted LD scores (Pan-UKB AMR), treat the results as
exploratory, and note that LDAK publishes no AMR tagging file (SumHer needs one
computed with `ldak --calc-tagging`). gencorr warns at start-up.

References per ancestry (all GRCh37 positions or rsID-keyed, which is what
gencorr matches on; check builds after download):

| ancestry | LDSC (`--ldsc_ld_<x>`)                                                                     | LAVA (`--lava_ref_<x>`)                                                | SumHer (`--sumher_tagfile_<x>`) |
| -------- | ------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------- | ------------------------------- |
| EUR      | `eur_w_ld_chr/` (Broad 1000G) or Pan-UKB EUR                                               | `lava-ukb-v1.1` (recommended by the LAVA authors) or 1000G `g1000_eur` | `ldak.thin.hapmap.gbr.tagging`  |
| EAS      | `eas_ldscores/` (Broad 1000G) or Pan-UKB EAS                                               | 1000G `g1000_eas` (LAVA `REFERENCE.md`)                                | `ldak.thin.hapmap.eas.tagging`  |
| AFR      | Pan-UKB AFR (covariate-adjusted) or 1000G AFR (Institut Pasteur, Zenodo 8096589)           | 1000G `g1000_afr`                                                      | `ldak.thin.hapmap.afr.tagging`  |
| SAS      | Pan-UKB CSA or 1000G SAS (Institut Pasteur)                                                | 1000G `g1000_sas`                                                      | `ldak.thin.hapmap.sas.tagging`  |
| AMR      | Pan-UKB AMR (covariate-adjusted; standard 1000G AMR scores are biased for admixed samples) | 1000G `g1000_amr` (exploratory)                                        | none published                  |

Pan-UKB LD scores (`UKBB.ALL.ldscore.tar.gz` from pan.ukbb.broadinstitute.org)
are one whole-genome file set per population, not split by chromosome: set
e.g. `--ldsc_ld_prefix_afr UKBB.AFR --ldsc_ld_split_by_chr_afr false`. The
`--ldsc_ld_prefix_<x>` / `--ldsc_ld_split_by_chr_<x>` options default to the
EUR ones. Use the same file set for `--ldsc_ld_<x>` and `--ldsc_w_ld_<x>`
(the default). 1000G PLINK references for LAVA (EUR, EAS, AFR, SAS, AMR) are
linked from LAVA's `REFERENCE.md`; the LAVA authors warn that small references
can inflate type I error, so treat non-EUR local results with care.

### Genome build and rsIDs

Every method matches variants **by rsID** and takes positions from its own
reference (LDSC LD scores, LAVA LD reference). The positions in your file are
never used for the analysis, so:

- **A file with rsIDs works on any build**, GRCh37 or GRCh38, without liftover.
  GWAS Catalog harmonised files (`hm_rsid`) are recognised automatically.
- **A file without rsIDs** (SNP IDs like `1:12345:A:G`, `chr1_12345_G_A`, or no
  SNP column at all) needs an **rsID map** for its build. The harmonizer then
  gives each such variant the rsID at the same chromosome and position with the
  same two alleles (either order, or the complementary strand). Variants that
  find no match are dropped and counted in `harmonized/<id>.qc.json`
  (`rsid_mapping`, `removed.no_rsid_after_mapping`).

|                     | GRCh37 files                                                                                         | GRCh38 files                                                                       |
| ------------------- | ---------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| samplesheet `build` | empty or `GRCh37`                                                                                    | `GRCh38`                                                                           |
| map parameter       | `--rsid_map_grch37`                                                                                  | `--rsid_map_grch38`                                                                |
| ready-made map      | the 1000 Genomes PLINK reference you use for LAVA: `--rsid_map_grch37 /refs/g1000_eur/g1000_eur.bim` | a GRCh38 `.bim` (e.g. 1000 Genomes high-coverage), or a table extracted from dbSNP |

A map is a PLINK `.bim`, or a table (tab/comma/space separated, optionally
gzipped) with a header containing CHR, BP and SNP columns and optionally A1/A2.
The usual names are recognised, so `#CHROM POS ID REF ALT` works, e.g. from a
VCF with rsIDs in the ID column (split multi-allelic sites first so each row
has one ALT allele):

```bash
bcftools norm -m -snps variants_GRCh38.vcf.gz -Ou \
  | bcftools query -i 'TYPE="snp" && ID~"^rs"' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
  | awk 'BEGIN{print "#CHROM\tPOS\tID\tREF\tALT"} {print}' | gzip > rsid_map_grch38.tsv.gz
```

Chromosomes must be named like in your summary statistics (`1`…`22` or
`chr1`…; a `chr` prefix is ignored). NCBI's dbSNP VCFs use RefSeq accessions
(`NC_000001.11`), so rename those first (`bcftools annotate --rename-chrs`).

Without alleles, matching uses chromosome + position only and skips positions
that carry several rsIDs. Only rows without an rsID are looked up; existing
rsIDs are never changed. Without a map, the QC report and the task log warn
when fewer than half of a file's IDs are rsIDs.

### GWAS Catalog studies

Put a study's accession in `sumstats` and gencorr downloads its **harmonised**
summary statistics from the GWAS Catalog FTP:

```csv
id,sumstats,n,ancestry
migraine,GCST90000016,550000,EUR
t2d_early,GCST000028,,EUR
mdd,/data/mdd.txt.gz,500199,EUR
```

- The file is found in `<range>/<accession>/harmonised/` (or, failing that, by
  an exact match in the Catalog's `harmonised_list.txt`). When a study has two
  harmonised files, the newer `<accession>.h.tsv.gz` is used.
- It is checked against the md5 in the study's metadata (`-meta.yaml`); the
  accession, URL, md5, metadata date, genome assembly, trait and samples are
  recorded under `source` in `harmonized/<id>.qc.json`.
- **N**: harmonised files rarely have an N column. gencorr uses, in order, the
  file's N column, the samplesheet `n` (or `n_cases + n_controls`), and the
  study's sample size from the metadata. The metadata size is used **only**
  when the study lists a single sample group, of the samplesheet's ancestry
  (European for `EUR`, East Asian for `EAS`). Otherwise the run stops and asks
  for `n`: a total over several ancestries does not match the LD references,
  and several groups can be stages or repeated entries (some studies list the
  same sample twice), so their sum is not reliably the analysed N. A study
  whose ancestry does not include the samplesheet's gets a warning in the task
  log and the QC report.
- Harmonised files carry rsIDs, so `build` and the rsID maps are not needed.
  Older files keep the submitted columns next to the harmonised `hm_` ones;
  gencorr then uses the `hm_` alleles and effects together (`colmap` overrides).
- Only harmonised studies can be given by accession. For others, download the
  file yourself and give its path or URL.
- Downloads need network access from the compute nodes. They are kept in the
  work directory (reused with `-resume`), or across runs with
  `--gwascatalog_cache <folder>` (one subfolder per accession).
- `--gwascatalog_base` points at a mirror or a local copy with the same layout.

### Column detection

The harmonizer recognises the usual names case-insensitively (full list in
`bin/harmonize_sumstats.py`, `COLUMN_SYNONYMS`):

| standard         | recognised source names (first match wins)                                       |
| ---------------- | -------------------------------------------------------------------------------- |
| SNP              | snp, rsid, rs_id, hm_rsid, variant_id, markername, snpid, snpid_ukb, …           |
| A1 (effect)      | a1, effect_allele, hm_effect_allele, ea, alt, allele1, tested_allele             |
| A2               | a2, other_allele, hm_other_allele, nea, ref, allele2, allele0, non_effect_allele |
| Z                | z, zscore, z_score, z-score, z_stat                                              |
| BETA             | beta, hm_beta, b, effect, logor, log_or, log_odds                                |
| OR               | or, odds_ratio, hm_odds_ratio (converted to BETA = ln OR)                        |
| P                | p, pval, p_value, p.value, …, p_bolt_lmm; else -log10 P from log10p / log10_p    |
| SE, N, EAF, INFO | the usual spellings (`r2` is not taken as INFO: it is often a model or LD r²)    |

`STAT` is **never** mapped automatically. Depending on the software it is a Z,
t or chi-square statistic, so map it with `colmap` if you know it is a Z.

**A1 must be the effect allele.** These layouts are recognised from their headers:

| software                              | effect allele      | other allele                                              |
| ------------------------------------- | ------------------ | --------------------------------------------------------- |
| GWAS Catalog harmonised               | `hm_effect_allele` | `hm_other_allele`                                         |
| METAL                                 | `Allele1`          | `Allele2`                                                 |
| SAIGE (has `AF_Allele2`/`AC_Allele2`) | `Allele2`          | `Allele1` (EAF from `AF_Allele2`)                         |
| PLINK 2 `--glm` (`A1`, `REF`, `ALT`)  | `A1`               | `OMITTED` if present, else whichever of REF/ALT is not A1 |
| REGENIE, BOLT-LMM                     | `ALLELE1`          | `ALLELE0`                                                 |

Any other file with only `Allele1`/`Allele2` columns is read the METAL way and
the QC report (and task log) carries `warning_allele_orientation`: check which
allele the effect refers to. Naming A1 or A2 in `colmap` turns the layout
detection off. Also check files where "A1" means something else (e.g. PLINK 1.9
`.assoc` A1 = minor allele, which _is_ the tested allele).

### Signed Z-score

Per variant, the first available of (`--z_priority`, default `z,p,beta_se`):

1. an explicit Z column;
2. |Z| from the two-sided P, sign from BETA/ln OR;
3. BETA / SE.

Option 2 before 3 matches what LDSC `munge_sumstats.py` and LAVA do with raw
files. That keeps results comparable.

## 2. Which pairs (`--pair_mode`)

| mode                       | pairs                                             |
| -------------------------- | ------------------------------------------------- |
| `all` (default)            | every pair of samplesheet traits                  |
| `target --targets mdd,scz` | each target vs every other trait (one vs group)   |
| `pairs --pairs pairs.csv`  | exactly the listed pairs (`trait1,trait2` header) |

Every method analyses the same pairs.

## 3. Reference data

### LDSC

| file                               | where                                                                          |
| ---------------------------------- | ------------------------------------------------------------------------------ |
| European LD scores `eur_w_ld_chr/` | Zenodo record 8182036 (`eur_w_ld_chr.tar.gz`, copy of the Alkes group archive) |
| HapMap3 SNP list `w_hm3.snplist`   | Zenodo record 7773502 (`w_hm3.snplist.gz`)                                     |

```bash
--ldsc_ld /refs/eur_w_ld_chr          # files 1.l2.ldscore.gz ... 22.l2.ldscore.gz
--ldsc_merge_alleles /refs/w_hm3.snplist
# use other ancestries' LD scores for non-European GWAS (match the GWAS ancestry!)
```

`--ldsc_merge_alleles` is strongly recommended. Without it, LDSC fails a whole
pair as soon as one rsID has different alleles in the two traits.

### LAVA

| reference                                                           | notes                                                                                                                             |
| ------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| **UK Biobank, LAVA LD format v1.1** (recommended for European data) | links in the LAVA repository's `REFERENCE.md`; ~15 GB unzipped, per chromosome; requires LAVA ≥ 0.1.5                             |
| 1000 Genomes PLINK files                                            | smaller and quicker, but the LAVA authors report possible type I error inflation and biased local h2 with 1000G for European data |

```bash
--lava_ref /refs/lava-ukb/lava-ukb-v1.1    # prefix WITHOUT the _chrN part
--lava_ref /refs/g1000/g1000_eur           # PLINK prefix (.bed/.bim/.fam)
```

All files starting with the prefix are staged, so keep the reference in its own
directory.

### SumHer (`--methods sumher`)

SumHer (LDAK 6.1) estimates global rg with a different heritability model
than LDSC. It needs an LDAK **tagging file** per ancestry:

| file                                                                       | where                                                                                                                                         |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `ldak.thin.hapmap.gbr.tagging` (EUR), `ldak.thin.hapmap.eas.tagging` (EAS) | dougspeed.com → Pre-computed Tagging Files (LDAK-Thin model, recommended by LDAK for genetic correlation; gunzip not needed, gencorr does it) |

```bash
--methods ldsc,lava,sumher \
--sumher_tagfile /refs/ldak/ldak.thin.hapmap.gbr.tagging.gz
```

SNP names in your data must match those in the tagging file. Check them with
`zcat ldak.thin.hapmap.gbr.tagging.gz | head -3`:

- **rsIDs** (`rs3094315`): the current pre-computed files (checked November
  2026). Nothing else to set.
- **`CHR:BP`** (`1:752566`): older releases and
  tagging files you made from a `CHR:BP`-named reference. Add
  `--sumher_snp_positions` with a GRCh37 rsID → position table (the 1000
  Genomes `.bim` you use for LAVA works, or a `CHR BP SNP` table), and gencorr
  renames each variant to `CHR:BP`.

With the wrong choice no SNP matches, and the pair gets `status = sumher_error`
with LDAK's message "contains summary statistics for none of the … predictors"
in `extra.error`. Options: `--check-sums NO --allow-ambiguous YES`, `--sumher_cutoff 0.01`,
`--sumher_genomic_control false` (LDAK 6's own default models multiplicative
inflation instead; set `true` for that).

LDAK reports rg and SEs but no p-value: gencorr gives z = rg / SE and a
two-sided normal p, marked in `extra.p_method`.

### Popcorn (`--methods popcorn`): cross-ancestry pairs

Popcorn (Brown et al. 2016) estimates the genetic correlation between a trait
measured in one population and a trait measured in another, so it runs **only on
cross-ancestry pairs** (EUR × EAS, EUR × AFR, EAS × SAS, ...); the other methods
run on the same-ancestry pairs (and Popcorn reports those as `same_ancestry`).
Popcorn assumes non-admixed samples, so pairs with an `AMR` trait get status
`admixed_ancestry`.

Each **ancestry pair needs a score file** (LD scores within and across the two
populations). gencorr takes the first of:

1. a row of `--popcorn_cfiles`, a CSV of pre-computed files:
   ```csv
   ancestry_1,ancestry_2,cfile
   EUR,EAS,/refs/popcorn/EUR_EAS.cscore
   EUR,AFR,/refs/popcorn/EUR_AFR.cscore
   ```
2. `--popcorn_cfile` with `--popcorn_cfile_ancestries` (default `EUR,EAS`): one
   pre-computed file, e.g. the Popcorn authors' 1000 Genomes EUR/EAS scores
   (link in Popcorn's README);
3. otherwise the score file is **computed** (`POPCORN_COMPUTE`, `popcorn
compute`) from one PLINK reference per ancestry: `--popcorn_ref` (EUR),
   `--popcorn_ref_eas`, `--popcorn_ref_afr`, `--popcorn_ref_sas`, e.g. the 1000
   Genomes `g1000_<x>` panels also used for LAVA. Restricted to
   `--popcorn_extract` (default: `--ldsc_merge_alleles`, i.e. HapMap3), with
   `--popcorn_maf` if set. Computing takes minutes for HapMap3 SNPs and longer
   for dense panels; the files are published to `popcorn/scores/` and kept
   across runs with `--popcorn_cache <folder>`. Population 1 is the ancestry
   that comes first in EUR, EAS, AFR, SAS.

   Popcorn's score files are **not bit-reproducible across machines**: the
   cross-population LD score of a few percent of SNPs moves by small steps
   with the CPU and the number of BLAS threads (on Popcorn's test panels: up to
   0.05 for 8% of SNPs; rg unchanged at 4 decimals). Compute each score file
   once, keep it (`--popcorn_cache`, or pass it with `--popcorn_cfiles`), and
   use the same file for all analyses you compare.

A pair with none of these stops the run at start-up, naming what is missing.
The score file's ancestry order is its population order (`--bfile1`,
`--bfile2`); gencorr puts each trait pair in that order.

```bash
--methods ldsc,lava,popcorn --popcorn_ref /refs/g1000_eur --popcorn_ref_afr /refs/g1000_afr \
    --popcorn_cache /refs/popcorn
```

- `--popcorn_gen_effect`: genetic effect instead of genetic impact correlation.
  The Popcorn authors call the effect model the more realistic one; the score
  file must be computed with `--gen_effect` too (computed files are, and their
  name records the model: `EUR_AFR.effect.cscore`).
- `--popcorn_use_mle`: maximum likelihood instead of the regression fit the
  authors recommend.
- Binary traits with `n_cases`, `n_controls` and `pop_prev` also get Popcorn's
  liability-scale h2 (`extra.h2_1_liab`, `h2_2_liab`).

**p-values:** Popcorn's own p for the correlation tests **rg < 1** (is the
cross-ancestry correlation lower than perfect?). gencorr's `p` column instead
tests **rg ≠ 0**, like every other method (z = rg / SE, normal), and keeps
Popcorn's test as `extra.z_pg_lt_1` / `extra.p_pg_lt_1`.

### LAVA loci

Loci: the default `--lava_loci` is LAVA's 2,495 GRCh37 blocks (fetched from the
pinned LAVA commit). **The locus file and the LAVA reference must be on the same
build** (both GRCh37 by default), because LAVA assigns SNPs to loci by the
reference's positions. Your summary statistics can be on any build as long as
they have rsIDs (see [Genome build and rsIDs](#genome-build-and-rsids)).

## 4. Important parameters

| parameter                                                                                                                                                         | default                                 |                                                                                                                                                                                                                                            |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--methods`                                                                                                                                                       | `ldsc,lava`                             | comma-separated: `ldsc`, `lava`, `sumher`, `popcorn`                                                                                                                                                                                       |
| `--lava_overlap`                                                                                                                                                  | `ldsc`                                  | sample-overlap correction from LDSC intercepts (LAVA vignette recipe); `none`; or a matrix file. The LAVA task stops if the part of the matrix it uses is not positive definite (e.g. \|r\| ≥ 1); LDSC's summary step warns about it first |
| `--lava_snp_set`                                                                                                                                                  | `pairwise`                              | `pairwise`: a pair uses SNPs shared by its two traits. `joint`: all traits loaded together, SNPs shared by **all** traits (LAVA's multi-trait behaviour)                                                                                   |
| `--lava_univ_thresh`                                                                                                                                              | `0.05`                                  | bivariate test only if both local h2 p-values are below this. For genome-wide screens the LAVA authors suggest a stricter (e.g. Bonferroni) filter                                                                                         |
| `--lava_loci_per_chunk`                                                                                                                                           | `250`                                   | parallelism only; results do not depend on it                                                                                                                                                                                              |
| `--lava_seed`                                                                                                                                                     | `20201231`                              | base seed for LAVA's Monte Carlo CIs/p-values                                                                                                                                                                                              |
| `--ldsc_rg_batch_size`                                                                                                                                            | `50`                                    | traits per `ldsc --rg` call; does not change estimates                                                                                                                                                                                     |
| `--drop_ambiguous`                                                                                                                                                | `false`                                 | remove A/T, C/G SNPs for all methods (LDSC munge always removes them)                                                                                                                                                                      |
| `--min_info`                                                                                                                                                      | `0` (off)                               | remove variants with INFO ≤ value, for all methods (see 4a)                                                                                                                                                                                |
| `--min_maf`                                                                                                                                                       | `0` (off)                               | remove variants with MAF ≤ value, for all methods (see 4a)                                                                                                                                                                                 |
| `--lava_snp_list`                                                                                                                                                 | none (all SNPs)                         | run LAVA only on the SNP IDs in this list, e.g. `w_hm3.snplist` (see 4b)                                                                                                                                                                   |
| `--sumher_tagfile`, `--sumher_tagfile_eas`                                                                                                                        | none                                    | LDAK tagging files for SumHer (EUR, EAS)                                                                                                                                                                                                   |
| `--sumher_snp_positions`                                                                                                                                          | none                                    | GRCh37 rsID → position table; SumHer variants named `CHR:BP`. Only for tagging files with `CHR:BP` names (older LDAK releases); current pre-computed files use rsIDs                                                                       |
| `--sumher_cutoff`, `--sumher_genomic_control`                                                                                                                     | `0.01`, `false`                         | LDAK `--cutoff`, `--genomic-control`                                                                                                                                                                                                       |
| `--popcorn_cfile`                                                                                                                                                 | none                                    | one Popcorn score file (default for EUR × EAS pairs)                                                                                                                                                                                       |
| `--popcorn_cfile_ancestries`                                                                                                                                      | `EUR,EAS`                               | populations of the score file, in order                                                                                                                                                                                                    |
| `--popcorn_cfiles`                                                                                                                                                | none                                    | CSV `ancestry_1,ancestry_2,cfile`: pre-computed score files for several ancestry pairs                                                                                                                                                     |
| `--popcorn_ref`, `--popcorn_ref_eas`, `--popcorn_ref_afr`, `--popcorn_ref_sas`                                                                                    | none                                    | PLINK references to compute score files for ancestry pairs without one                                                                                                                                                                     |
| `--popcorn_extract`, `--popcorn_maf`, `--popcorn_cache`                                                                                                           | `--ldsc_merge_alleles`, Popcorn's, none | SNP list, MAF cut-off and cache folder for computed score files                                                                                                                                                                            |
| `--popcorn_gen_effect`, `--popcorn_use_mle`                                                                                                                       | `false`, `false`                        | genetic effect correlation; MLE fit                                                                                                                                                                                                        |
| `--ldsc_ld_eas`, `--ldsc_w_ld_eas`, `--lava_ref_eas`                                                                                                              | none                                    | references for traits with `ancestry = EAS` (see [Ancestry](#ancestry))                                                                                                                                                                    |
| `--ldsc_ld_<x>`, `--ldsc_w_ld_<x>`, `--lava_ref_<x>`, `--sumher_tagfile_<x>` (x = afr, sas, amr); `--ldsc_ld_prefix_<x>`, `--ldsc_ld_split_by_chr_<x>` (also eas) | none                                    | references for traits with `ancestry = AFR/SAS/AMR`; prefix/layout default to the EUR settings (see [Ancestry](#ancestry))                                                                                                                 |
| `--rsid_map_grch37`, `--rsid_map_grch38`                                                                                                                          | none                                    | rsID maps for files without rsIDs, per samplesheet `build` (see [Genome build and rsIDs](#genome-build-and-rsids))                                                                                                                         |
| `--gwascatalog_cache`                                                                                                                                             | none                                    | folder that keeps GWAS Catalog downloads across runs (see [GWAS Catalog studies](#gwas-catalog-studies))                                                                                                                                   |
| `--gwascatalog_base`                                                                                                                                              | GWAS Catalog FTP                        | mirror or local copy of the GWAS Catalog summary statistics                                                                                                                                                                                |
| `--save_harmonized` / `--save_intermediates`                                                                                                                      | `false`                                 | publish intermediate files                                                                                                                                                                                                                 |

`nextflow run mauricioroza/gencorr --help --show_hidden` lists everything.

## 4a. Variant QC: `--min_info`, `--min_maf`

Applied once, during harmonization, so **every method** sees the same
filtered variants. The rule is LDSC munge's: a variant is removed when its value
is **≤** the threshold. Missing values are kept (genotyped SNPs often have no
INFO).

- INFO comes from the INFO column (`info`, `info_ukb`, `impinfo`, `r2`, …).
- MAF comes from a `MAF` column if present, otherwise `min(EAF, 1 − EAF)`.
- If a file has no such column the threshold cannot be applied: nothing is
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
(first column; header optional; plain or gzipped).

|                                  | all SNPs (default)                           | HapMap3 only (`--lava_snp_list w_hm3.snplist`)       |
| -------------------------------- | -------------------------------------------- | ---------------------------------------------------- |
| information per locus            | all SNPs; LAVA compresses LD to PCs itself   | ~15–20% of SNPs; some loci drop to a handful of SNPs |
| robustness to imputation quality | depends on QC → use `--min_info`/`--min_maf` | common, well-imputed SNPs only                       |
| same SNPs across GWAS            | pairwise intersection                        | fixed list                                           |
| ancestry                         | whatever the reference covers                | HapMap3 is less complete outside European ancestry   |

LDSC is not affected (it has its own list: `--ldsc_merge_alleles`). The SNP set
used is recorded in every LAVA row (`extra.snp_set`: `all` or the file name).

## 5. Running on a SLURM cluster

```bash
# Nextflow >= 25.10.4 and Apptainer (or Singularity) on the login node
export NXF_HOME=/path/on/project/storage/$USER/.nextflow
export NXF_APPTAINER_CACHEDIR=/path/on/project/storage/$USER/apptainer_cache

nextflow run mauricioroza/gencorr -profile apptainer -c cluster.config \
    --input samplesheet.csv --outdir results ...
```

with a small `cluster.config` (adjust account, partition and limits to your
allocation):

```groovy
process {
    executor       = 'slurm'
    clusterOptions = '-A <project_account>'
    resourceLimits = [ cpus: 16, memory: 120.GB, time: 240.h ]
}
executor { queueSize = 200; submitRateLimit = '10/1min' }
```

Many clusters have a ready-made institutional profile in
[nf-core/configs](https://github.com/nf-core/configs) (`-profile <name>`,
loaded automatically from `--custom_config_base`); check that it matches your
system before relying on it.

The `gencorr-tools` image is private (LAVA's licence): on a cluster without
Docker, build it elsewhere and convert it once,
`apptainer build gencorr-tools_0.1.0.sif docker-daemon://ghcr.io/mauricioroza/gencorr-tools:0.1.0`,
or pull it from ghcr.io with a token that has access. Compute nodes need
network access for GWAS Catalog accessions and the default `--lava_loci` URL
(or give local copies).

Memory: harmonization loads one GWAS file at a time in pandas (`process_low`,
12 GB; a 9.5 M-row file peaks at about 4.4 GB). For very large files
(> 20 M rows) raise it in a config:
`process { withName: 'HARMONIZE_SUMSTATS' { memory = 32.GB } }`.

`-resume` works as usual. Changing a LAVA parameter re-runs only LAVA.

## 6. Running on a cloud VM

[cloud_oci.md](cloud_oci.md) sets up one Oracle Cloud VM that runs the pipeline
with the local executor and Docker (scripts in `cloud/oci/`), with the reasons
for each choice.
