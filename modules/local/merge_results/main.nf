/*
 * MERGE_RESULTS
 * -------------
 * Concatenates standard results tables, validates them against the gencorr
 * contract (bin/gencorr_common.py: RESULT_COLUMNS), rejects duplicated
 * estimates and sorts. Used for LAVA chunks and for the final all-method table.
 */
process MERGE_RESULTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/karlssonlab/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(tables, stageAs: 'in/table_*.tsv')

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: results
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_results.py --out ${prefix}.tsv ${tables}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
