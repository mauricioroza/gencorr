# Comparing gencorr with BigaGWAS on the test dataset

Goal: run the **same files** through BigaGWAS (bigagwas.org) and gencorr, with
the **same references and settings**, and check that the numbers agree.

## 0. Why the `-profile test` output cannot be compared directly

The test profile uses small test references (LD scores computed from 42
regions, a 1000G subset). BigaGWAS uses the full references:

| | BigaGWAS (from its source code) | gencorr `-profile test` |
|---|---|---|
| LDSC LD scores | `eur_w_ld_chr` | `tests/data/ldsc` (test only) |
| LDSC SNP list | `w_hm3.snplist` | test allele list |
| LAVA reference | 1000 Genomes EUR (PLINK) | `g1000_test` (vignette subset) |
| LAVA sample overlap | one number per job, default 0 | from LDSC intercepts |
| LAVA case/control | not asked in the 2-trait form | from `n_cases`/`n_controls` |

So gencorr is re-run below with BigaGWAS' references and settings.

## 1. Files to upload

`tests/bigagwas/prepare_upload.py` converts the 4 vignette files into what
BigaGWAS accepts (it needs `a1 a2 se pvalue n` + `beta` or `OR`; the raw
`neuro`, `depression` and `rheuma` files would be rejected). The same files are
used as gencorr input, so both tools start from identical data.

| file | trait | effect column | notes |
|---|---|---|---|
| `bmi_biga.txt.gz` | Body mass index | BETA | |
| `neuro_biga.txt.gz` | Neuroticism | BETA | Z column removed: both tools derive Z from P |
| `depression_biga.txt.gz` | Major depression | OR | SE = SE of log OR |
| `rheuma_biga.txt.gz` | Rheumatoid arthritis | OR | SE derived from OR and P (raw file has none); OR = 1.00 variants kept with SE at the rounding bound; 391 variants with P = 1 dropped |

Regenerate with:
`python3 tests/bigagwas/prepare_upload.py <LAVA vignettes/data dir> biga_upload/`

## 2. BigaGWAS jobs (12)

Use the **2 traits** analysis. One job per pair and method:

| pairs (6) | bmi–neuro, bmi–depression, bmi–rheuma, neuro–depression, neuro–rheuma, depression–rheuma |
|---|---|

For every job:
* Upload both files; genome build **hg19 / GRCh37**; population **European**.
* Column names should be detected automatically (SNP, CHR, BP, A1, A2, BETA or OR, SE, P, N).
  Check the job log's "Interpreting columns as" section.
* Name the traits exactly `bmi`, `neuro`, `depression`, `rheuma`. Then no
  trait map is needed when comparing.

**LDSC jobs:** method LDSC; leave the liability (sample/population prevalence)
fields empty. rg does not depend on them, and the comparison uses rg.

**LAVA jobs:** method LAVA with
`univ_thresh = 0.05`, p-values **yes**, confidence intervals **yes**,
`param_lim = 1.25`, sample overlap **0**.

Write down each job ID; keep the result downloads.

## 3. gencorr with BigaGWAS' references

Download once (WSL or cluster):

```bash
mkdir -p ~/refs && cd ~/refs
wget -O eur_w_ld_chr.tar.gz 'https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1'
tar xzf eur_w_ld_chr.tar.gz                       # -> eur_w_ld_chr/1..22.l2.ldscore.gz
wget -O w_hm3.snplist.gz 'https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1'
gunzip w_hm3.snplist.gz
wget -O g1000_eur.zip https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download
unzip g1000_eur.zip                               # check the prefix: ls g1000_eur*
```

Check the extracted names (`ls eur_w_ld_chr | head`, `ls g1000_eur*`) and
adjust the paths below if they differ.

```bash
cd ~/gencorr
nextflow run . -profile docker \
    --input tests/bigagwas/samplesheet_biga.csv \
    --outdir biga_compare \
    --methods ldsc,lava \
    --ldsc_ld ~/refs/eur_w_ld_chr --ldsc_merge_alleles ~/refs/w_hm3.snplist \
    --lava_ref ~/refs/g1000_eur/g1000_eur \
    --lava_loci https://raw.githubusercontent.com/josefin-werme/LAVA/e729a245f7b6923967a96804fbf5246eadf2d6c6/vignettes/data/test.loci \
    --lava_overlap none \
    --lava_snp_set pairwise
```

* `samplesheet_biga.csv` points to `biga_upload/…`: put the `biga_upload/`
  folder inside `~/gencorr/` (paths are relative to where you launch).
* It has no `n_cases`/`n_controls`, so LAVA treats every trait as continuous,
  as BigaGWAS' 2-trait form appears to (it does not ask for cases/controls).
* `--lava_loci` = the 42 vignette loci. All of them are standard LAVA blocks,
  so their IDs match BigaGWAS' 2,495-locus output. Loci outside these 42 have no
  data in these files anyway.

## 4. Enter BigaGWAS results and compare

**LDSC**: fill `tests/bigagwas/biga_ldsc_template.csv` (one row per job; copy
rg, se, z, p from the job's result log, "Summary of Genetic Correlation
Results" table), save as `biga_ldsc.csv`, then:

```bash
bin/compare_results.py --test biga_compare/gencorr_results.tsv --method ldsc \
    --reference biga_ldsc.csv --reference-format custom \
    --colmap trait1=trait1,trait2=trait2,rg=rg,rg_se=se,z=z,p=p \
    --report ldsc_vs_biga
```

**LAVA**: combine the 6 BigaGWAS LAVA result tables into one file with a
column for each trait name (`trait1`, `trait2`) and the locus ID, then map the
columns BigaGWAS uses (look at the header first; LAVA's own names are `locus`,
`phen1`, `phen2`, `rho`, `rho.lower`, `rho.upper`, `p`):

BigaGWAS' `<job>_out_bivar.txt` is comma-separated and uses the trait names
typed in the job; convert and map them first, e.g.:

```bash
sed 's/,/\t/g' 971_out_bivar.txt > biga_lava.tsv
printf 'reference_name,gencorr_id\ntype_II_diabetes,t2d\ncoronary_heart_disease,chd\n' > names.csv
# then add --trait-map names.csv to the command below
```

```bash
bin/compare_results.py --test biga_compare/gencorr_results.tsv --method lava --only-ok \
    --reference biga_lava.tsv --reference-format custom \
    --colmap trait1=phen1,trait2=phen2,locus=locus,rg=rho,rg_ci_low=rho.lower,rg_ci_high=rho.upper,p=p \
    --tol rg=1e-3,p=0.2,rg_ci_low=0.05,rg_ci_high=0.05 \
    --allow-missing --report lava_vs_biga
```

`ldsc_vs_biga.txt` / `lava_vs_biga.txt` summarise the agreement;
`*.pairs.tsv` lists every pair with both values side by side.

### LAVA: use BigaGWAS' SNP set

BigaGWAS runs LAVA on a HapMap3-based SNP list only. Add
`--lava_snp_list ~/refs/w_hm3.snplist` to the gencorr LAVA run; without it the
inputs differ (~5× more SNPs per locus) and so do the estimates.

With the list, gencorr reproduces BigaGWAS' LAVA results (section 6).
BigaGWAS uses LAVA 0.1.0 and gencorr 0.1.5; for continuous traits on the same
SNPs this made no visible difference.

## 5. What agreement to expect

| quantity | expectation | if it differs |
|---|---|---|
| LDSC rg, se, z, p | equal to the 4th decimal (same LDSC 1.0.1, same references) | BigaGWAS harmonizes against its own reference VCF first (drops or flips some variants). Compare the "SNPs with valid alleles" counts in both logs |
| LAVA rho | close (≤ 1e-3) | different SNP set after BigaGWAS harmonization; a different LAVA version on their server; case/control handling |
| LAVA p, CIs | similar but not identical | Monte Carlo; different seeds. Compare on the log scale |
| which loci are tested | same `univ_filtered` pattern | univariate p near the 0.05 threshold can flip with tiny input differences |

Differences that are not explained by this table are worth reporting, either
to us or to the BigaGWAS team.

## 6. Results on the test jobs

gencorr run with BigaGWAS' references and settings (sections 3–4, LAVA with
`--lava_snp_list w_hm3.snplist`), compared with `bin/compare_results.py`.
Jobs 971/972 were re-checked on 2026-11-08 against all 2,495 standard loci
(**PASS**, tolerances as in section 4); 966/967 are from the first comparison
(42 vignette loci).

| BigaGWAS job | traits | method | result |
|---|---|---|---|
| 966 | bmi × depression | LDSC | rg 0.1154 (gencorr) vs 0.1153; BigaGWAS' harmonization used 4 fewer SNPs |
| 967 | bmi × depression | LAVA | 32/32 loci in common; SNPs and PCs identical in 29/32 (3 differ by 1–2 SNPs BigaGWAS dropped); rho median difference 1e-5; same loci at p < 0.05 |
| 971 | T2D × CHD | LAVA | same 8 loci tested; SNPs and PCs identical 8/8; rho within 1e-7; same 2 loci at p < 0.05 (1126, 1801); log10 p within 0.005, CI bounds within 0.019 (Monte Carlo); univariate h2 identical (63/63 locus × trait) |
| 972 | T2D × CHD | LDSC | identical: rg 0.3476, se 0.0796, z 4.3684, p 1.2518e-05, 2,937 SNPs; h2, intercepts and gencov identical in both logs |

Without the HapMap3 list (gencorr's default, all SNPs), LAVA rho on job 967
correlates r = 0.92 with BigaGWAS; the difference is entirely the SNP set.

The T2D (ieu-a-25) and CHD (ieu-a-9) example files are Metabochip
(regional) data, not genome-wide. They test that the computations agree, but
the LDSC estimates themselves are not meaningful (e.g. T2D h2 = 5.6).
