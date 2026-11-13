# mauricioroza/gencorr: Citations

## How to cite

gencorr has no DOI yet. Until its first release, cite the repository and the
version you ran (`pipeline_info/gencorr_software_versions.yml`, entry
`Workflow`), e.g. *mauricioroza/gencorr v0.1.0dev,
https://github.com/mauricioroza/gencorr*.

Please also cite the **methods** you ran and the **reference data** you used
(sections below), and the nf-core and Nextflow papers.

## Pipeline framework

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Methods

Cite the method papers for the methods in `--methods`.

### LDSC (`--methods ldsc`)

- [LD score regression: SNP heritability](https://doi.org/10.1038/ng.3211)

  > Bulik-Sullivan BK, Loh PR, Finucane HK, Ripke S, Yang J; Schizophrenia Working Group of the Psychiatric Genomics Consortium, Patterson N, Daly MJ, Price AL, Neale BM. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2015;47(3):291-295. doi: 10.1038/ng.3211.

- [Cross-trait LD score regression: genetic correlation](https://doi.org/10.1038/ng.3406)

  > Bulik-Sullivan B, Finucane HK, Anttila V, et al. An atlas of genetic correlations across human diseases and traits. Nat Genet. 2015;47(11):1236-1241. doi: 10.1038/ng.3406.

- [Liability-scale heritability](https://doi.org/10.1016/j.ajhg.2011.02.002), when `n_cases`, `n_controls` and `pop_prev` are given

  > Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. Am J Hum Genet. 2011;88(3):294-305. doi: 10.1016/j.ajhg.2011.02.002.

### LAVA (`--methods lava`)

- [LAVA: local genetic correlation](https://doi.org/10.1038/s41588-022-01017-y). Also the source of the default locus definitions (`--lava_loci`, 2,495 GRCh37 blocks) and of the sample-overlap correction from LDSC intercepts (`--lava_overlap ldsc`).

  > Werme J, van der Sluis S, Posthuma D, de Leeuw CA. An integrated framework for local genetic correlation analysis. Nat Genet. 2022;54(3):274-282. doi: 10.1038/s41588-022-01017-y.

### SumHer (`--methods sumher`)

- [SumHer / LDAK](https://doi.org/10.1038/s41588-018-0279-5)

  > Speed D, Balding DJ. SumHer better estimates the SNP heritability of complex traits from summary statistics. Nat Genet. 2019;51(2):277-284. doi: 10.1038/s41588-018-0279-5.

  The LDAK pre-computed tagging files were built from UK Biobank data; also follow UK Biobank's acknowledgement policy when you use them.

### Popcorn (`--methods popcorn`)

- [Popcorn: trans-ethnic genetic correlation](https://doi.org/10.1016/j.ajhg.2016.05.001)

  > Brown BC; Asian Genetic Epidemiology Network Type 2 Diabetes Consortium, Ye CJ, Price AL, Zaitlen N. Transethnic genetic-correlation estimates from summary statistics. Am J Hum Genet. 2016;99(1):76-88. doi: 10.1016/j.ajhg.2016.05.001.

## Reference data

Cite the data behind the references you pass to the pipeline.

- [1000 Genomes Project](https://doi.org/10.1038/nature15393): the LDSC European LD scores (`eur_w_ld_chr`, `--ldsc_ld`) and the 1000 Genomes LAVA reference (`g1000_eur`, `--lava_ref`)

  > The 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature. 2015;526(7571):68-74. doi: 10.1038/nature15393.

- [HapMap 3](https://doi.org/10.1038/nature09298): the SNP list `w_hm3.snplist` (`--ldsc_merge_alleles`, `--lava_snp_list`)

  > The International HapMap 3 Consortium. Integrating common and rare genetic variation in diverse human populations. Nature. 2010;467(7311):52-58. doi: 10.1038/nature09298.

- [UK Biobank](https://doi.org/10.1038/s41586-018-0579-z): the LAVA UK Biobank LD reference (`lava-ukb-v1.1`, `--lava_ref`); also follow UK Biobank's own acknowledgement policy

  > Bycroft C, Freeman C, Petkova D, et al. The UK Biobank resource with deep phenotyping and genomic data. Nature. 2018;562(7726):203-209. doi: 10.1038/s41586-018-0579-z.

## Related work

- [BIGA](https://doi.org/10.1101/2023.04.28.538585) (bigagwas.org): cloud platform for cross-trait genetic analyses; its results were used to compare and test gencorr's results ([docs/bigagwas_comparison.md](docs/bigagwas_comparison.md)).

  > Li Y, Xue F, et al., Zhao B. Analyzing bivariate cross-trait genetic architecture in GWAS summary statistics with the BIGA cloud computing platform. bioRxiv. 2023. doi: 10.1101/2023.04.28.538585.

## Software used inside the pipeline

The harmonization, export and parsing scripts in `bin/` (Python) and the LAVA
runner (`bin/run_lava.R`):

- [pandas](https://doi.org/10.25080/Majora-92bf1922-00a)

  > McKinney W. Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference. 2010:56-61. doi: 10.25080/Majora-92bf1922-00a.

- [NumPy](https://doi.org/10.1038/s41586-020-2649-2)

  > Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585(7825):357-362. doi: 10.1038/s41586-020-2649-2.

- [SciPy](https://doi.org/10.1038/s41592-019-0686-2)

  > Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 2020;17(3):261-272. doi: 10.1038/s41592-019-0686-2.

- [R](https://www.R-project.org/)

  > R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

## Testing

- [nf-test](https://doi.org/10.1093/gigascience/giaf130): module, function and pipeline tests ([docs/validation.md](docs/validation.md))

  > Forer L, Schönherr S. Improving the reliability, quality, and maintainability of bioinformatics pipelines with nf-test. GigaScience. 2025. doi: 10.1093/gigascience/giaf130.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/): LDSC 1.0.1 package

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/): LDSC container image

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity / Apptainer](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
