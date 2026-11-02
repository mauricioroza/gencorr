/*
 * EXPORT_SUMSTATS
 * ---------------
 * Writes a harmonized file in the column layout one method expects
 * (bin/export_sumstats.py: FORMATS). Generic: the format is an input, so a new
 * method reuses this module by adding a FORMATS entry.
 *
 * Output name: <id>.<format>.tsv (or .tsv.gz for compressed formats)
 */
process EXPORT_SUMSTATS {
    tag "${meta.id}:${format}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/karlssonlab/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(harmonized)
    val format

    output:
    tuple val(meta), path("${meta.id}.${format}.tsv*"), emit: sumstats
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_sumstats.py \\
        --input ${harmonized} \\
        --format ${format} \\
        --out ${meta.id}.${format}.tsv
    """

    stub:
    """
    touch ${meta.id}.${format}.tsv
    """
}
