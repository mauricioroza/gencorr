/*
 * HARMONIZE_SUMSTATS
 * ------------------
 * Converts one raw GWAS summary statistics file into the gencorr harmonized
 * format (bin/gencorr_common.py: HARMONIZED_COLUMNS) and writes a QC report.
 * All method-specific inputs are derived from this one file, so every method
 * sees the same variants, alleles and Z-scores.
 *
 * meta fields used: id, n, n_cases, n_controls, colmap (all but id optional)
 * ext.args: extra harmonize_sumstats.py flags, e.g. '--drop-ambiguous'
 */
process HARMONIZE_SUMSTATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/genis/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(sumstats)

    output:
    tuple val(meta), path("${prefix}.harmonized.tsv.gz"), emit: sumstats
    tuple val(meta), path("${prefix}.qc.json")          , emit: qc
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'")          , emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def n_arg    = meta.n          ? "--n ${meta.n}"                   : ''
    def case_arg = meta.n_cases    ? "--n-cases ${meta.n_cases}"       : ''
    def ctrl_arg = meta.n_controls ? "--n-controls ${meta.n_controls}" : ''
    def colmap   = meta.colmap     ? "--colmap '${meta.colmap}'"       : ''
    """
    harmonize_sumstats.py \\
        --input ${sumstats} \\
        --id ${meta.id} \\
        --out ${prefix} \\
        ${n_arg} ${case_arg} ${ctrl_arg} ${colmap} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'SNP\\tCHR\\tBP\\tA1\\tA2\\tZ\\tBETA\\tSE\\tP\\tN\\tEAF\\tINFO\\n' | gzip > ${prefix}.harmonized.tsv.gz
    echo '{}' > ${prefix}.qc.json
    """
}
