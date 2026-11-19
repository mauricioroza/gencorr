/*
 * GWAMA_COMBINE
 * -------------
 * N-weighted multivariate GWAMA (Baselmans et al. 2019): combines the
 * harmonized files of a group of genetically correlated traits into one new
 * trait in the same harmonized format (bin/gwama.py).
 *
 * meta fields used: id (name of the combined trait), gwama_traits (trait IDs,
 *   in the order of the sumstats files), gwama_signs (1 or -1 per trait)
 * ldsc_h2, ldsc_results: LDSC_SUMMARISE tables for (at least) the group's
 *   traits and all their pairs: observed-scale h2, univariate intercepts and
 *   cross-trait intercepts (gcov_int)
 * ext.args: extra gwama.py flags, e.g. '--n-out observed --min-rg 0.5'
 */
process GWAMA_COMBINE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(sumstats, stageAs: 'input/*')
    path ldsc_h2
    path ldsc_results

    output:
    tuple val(meta), path("${prefix}.harmonized.tsv.gz"), emit: sumstats
    tuple val(meta), path("${prefix}.gwama.json")       , emit: qc
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'")          , emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('numpy'),  eval("python3 -c 'import numpy; print(numpy.__version__)'")            , emit: versions_numpy, topic: versions
    tuple val("${task.process}"), val('scipy'),  eval("python3 -c 'import scipy; print(scipy.__version__)'")            , emit: versions_scipy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def signs = meta.gwama_signs ? "--signs ${meta.gwama_signs.join(' ')}" : ''
    """
    gwama.py \\
        --id ${meta.id} \\
        --traits ${meta.gwama_traits.join(' ')} \\
        --sumstats ${sumstats.join(' ')} \\
        ${signs} \\
        --ldsc-h2 ${ldsc_h2} \\
        --ldsc-results ${ldsc_results} \\
        --out ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'SNP\\tCHR\\tBP\\tA1\\tA2\\tZ\\tBETA\\tSE\\tP\\tN\\tEAF\\tINFO\\n' | gzip > ${prefix}.harmonized.tsv.gz
    echo '{}' > ${prefix}.gwama.json
    """
}
