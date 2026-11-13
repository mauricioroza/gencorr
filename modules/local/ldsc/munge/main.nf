/*
 * LDSC_MUNGE
 * ----------
 * munge_sumstats.py on the LDSC export of a harmonized file. Applies LDSC's own
 * QC (INFO > 0.9, MAF > 0.01 when those columns are present, strand-ambiguous
 * and non-SNP removal) and, if given, restricts/aligns to --merge-alleles
 * (normally w_hm3.snplist).
 *
 * The Z column carries the sign; munge recomputes |Z| from P.
 *
 * --chunksize 500000: munge reads the file in chunks (default 5M rows) and
 * filters each with Series.isin(merge-alleles SNPs). In this container's
 * pandas (0.20.3), isin on more than 1M rows switches to np.in1d, which for
 * text compares every SNP with every merge-alleles SNP: a 9M-row file took
 * hours instead of a minute. Chunks of <= 1M rows keep the hash-table path.
 * The per-chunk filters are row-wise, so results are identical.
 */
process LDSC_MUNGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(sumstats)
    path merge_alleles   // optional: pass [] for none

    output:
    tuple val(meta), path("${meta.id}.sumstats.gz"), emit: sumstats
    tuple val(meta), path("${meta.id}.munge.log")  , emit: log
    tuple val("${task.process}"), val('ldsc'), eval("grep -m1 -o 'Version [0-9.]*' ${meta.id}.munge.log | cut -d' ' -f2"), emit: versions_ldsc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def merge = merge_alleles ? "--merge-alleles ${merge_alleles}" : ''
    """
    munge_sumstats.py \\
        --sumstats ${sumstats} \\
        --out ${meta.id} \\
        --signed-sumstats Z,0 \\
        --snp SNP --a1 A1 --a2 A2 --p P --N-col N \\
        ${merge} \\
        --chunksize 500000 \\
        ${args}
    mv ${meta.id}.log ${meta.id}.munge.log
    """

    stub:
    """
    echo "SNP A1 A2 Z N" | gzip > ${meta.id}.sumstats.gz
    printf '* Version 1.0.1\\n' > ${meta.id}.munge.log
    """
}
