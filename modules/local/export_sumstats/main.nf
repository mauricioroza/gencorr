/*
 * EXPORT_SUMSTATS
 * ---------------
 * Writes a harmonized file in the column layout one method expects
 * (bin/export_sumstats.py: FORMATS). Generic: the format is an input, so a new
 * method reuses this module by adding a FORMATS entry.
 *
 * Output name: <id>.<format>.tsv (or .tsv.gz for compressed formats)
 *
 * snp_list (optional, [] for none): keep only these SNP IDs (e.g. w_hm3.snplist).
 *   For format 'sumher' it is instead an rsID -> position table (PLINK .bim or
 *   CHR/BP/SNP table): variants are then named CHR:BP, for tagging files with
 *   CHR:BP names (older LDAK releases; current ones use rsIDs).
 */
process EXPORT_SUMSTATS {
    tag "${meta.id}:${format}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(harmonized)
    val format
    path snp_list    // optional: [] for none

    output:
    tuple val(meta), path("${meta.id}.${format}.tsv*"), emit: sumstats
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'")          , emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def restrict = !snp_list ? '' :
                   format == 'sumher' ? "--predictor-positions ${snp_list}" : "--snp-list ${snp_list}"
    """
    export_sumstats.py \\
        --input ${harmonized} \\
        --format ${format} \\
        --out ${meta.id}.${format}.tsv \\
        ${restrict} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}.${format}.tsv
    """
}
