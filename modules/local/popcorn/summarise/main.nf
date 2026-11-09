/*
 * POPCORN_SUMMARISE
 * -----------------
 * Parses all Popcorn outputs into the standard results table
 * (bin/parse_popcorn.py): one row per cross-ancestry pair, failed pairs kept
 * as 'popcorn_error'. The p column tests rg != 0; Popcorn's own rg < 1 test
 * is kept in extra (z_pg_lt_1, p_pg_lt_1).
 */
process POPCORN_SUMMARISE {
    tag "popcorn"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    path results, stageAs: 'results/*'
    path logs   , stageAs: 'logs/*'

    output:
    path "gencorr.popcorn.tsv", emit: results
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'")          , emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('scipy'),  eval("python3 -c 'import scipy; print(scipy.__version__)'")            , emit: versions_scipy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    parse_popcorn.py --results ${results} --logs ${logs} --out gencorr ${args}
    """

    stub:
    """
    touch gencorr.popcorn.tsv
    """
}
