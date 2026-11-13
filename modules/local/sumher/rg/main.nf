/*
 * SUMHER_RG
 * ---------
 * SumHer (LDAK) global genetic correlation for one trait pair:
 *     ldak6 --sum-cors <pair> --summary t1 --summary2 t2 --tagfile <tagging>
 * Options
 * (--check-sums, --allow-ambiguous, --cutoff, --genomic-control) come from
 * ext.args (conf/modules.config).
 *
 * A failing LDAK run does not stop the pipeline: the log is kept and an empty
 * .cors is written, and SUMHER_SUMMARISE reports the pair with status
 * 'sumher_error' and LDAK's message (as LDSC errors are reported).
 *
 * meta    : [ id: '<trait1>__<trait2>', trait1, trait2, ancestry ]
 * sumstats: SumHer-format files named <trait id>.sumher.tsv
 * tagfile : LDAK tagging file (plain or .gz) whose predictor names match the
 *           summary statistics (rsIDs, or CHR:BP for older LDAK tagging files)
 */
process SUMHER_RG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-ldak:6.1"

    input:
    tuple val(meta), path(sumstats1), path(sumstats2)
    path tagfile

    output:
    tuple val(meta), path("${meta.id}.cors"), emit: cors
    tuple val(meta), path("${meta.id}.log") , emit: log
    tuple val("${task.process}"), val('ldak'), eval("ls \$(dirname \$(command -v ldak6))/../conda-meta/ldak6-*.json | sed -E 's/.*ldak6-([0-9.]+)-.*/\\1/'"), emit: versions_ldak, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def gz     = tagfile.name.endsWith('.gz')
    def tag    = gz ? tagfile.name[0..-4] : tagfile.name
    """
    ${gz ? "gunzip -c ${tagfile} > ${tag}" : ''}
    if ! ldak6 --sum-cors ${meta.id} \\
            --summary ${sumstats1} \\
            --summary2 ${sumstats2} \\
            --tagfile ${tag} \\
            ${args} > ${meta.id}.log 2>&1; then
        echo "LDAK failed for ${meta.id}; reported as sumher_error" >&2
        tail -n 5 ${meta.id}.log >&2
    fi
    touch ${meta.id}.cors
    """

    stub:
    def args   = task.ext.args ?: ''
    // the LDAK options are recorded so -stub runs can check them (tests/python/test_cli_options.py)
    """
    printf -- '--summary ${sumstats1}\\n--summary2 ${sumstats2}\\n${args}\\n' > ${meta.id}.log
    touch ${meta.id}.cors
    """
}
