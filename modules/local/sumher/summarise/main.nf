/*
 * SUMHER_SUMMARISE
 * ----------------
 * Parses all SumHer (LDAK --sum-cors) outputs into the standard results table
 * (bin/parse_sumher.py): one row per pair, failed pairs kept as 'sumher_error'.
 */
process SUMHER_SUMMARISE {
    tag "sumher"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    path logs, stageAs: 'logs/*'
    path cors, stageAs: 'cors/*'

    output:
    path "gencorr.sumher.tsv", emit: results
    tuple val("${task.process}"), val('python'), eval("python3 -c 'import platform; print(platform.python_version())'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'")          , emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('scipy'),  eval("python3 -c 'import scipy; print(scipy.__version__)'")            , emit: versions_scipy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    parse_sumher.py --logs ${logs} --cors ${cors} --out gencorr ${args}
    """

    stub:
    """
    touch gencorr.sumher.tsv
    """
}
