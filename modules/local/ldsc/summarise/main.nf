/*
 * LDSC_SUMMARISE
 * --------------
 * Parses all LDSC h2 and rg logs into:
 *   ldsc.tsv              standard results (one row per pair)
 *   ldsc_h2.tsv           per-trait heritability, intercept, ratio, lambda GC
 *   sample_overlap.txt    sampling-correlation matrix for LAVA (LAVA vignette recipe)
 */
process LDSC_SUMMARISE {
    tag "ldsc"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/genis/gencorr-tools:0.1.0"

    input:
    path h2_logs, stageAs: 'h2/*'
    path rg_logs, stageAs: 'rg/*'

    output:
    path "gencorr.ldsc.tsv"          , emit: results
    path "gencorr.ldsc_h2.tsv"       , emit: h2
    path "gencorr.sample_overlap.txt", emit: overlap
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rg = rg_logs ? "--rg-logs ${rg_logs}" : ''
    """
    parse_ldsc.py --h2-logs ${h2_logs} ${rg} --out gencorr
    """

    stub:
    """
    touch gencorr.ldsc.tsv gencorr.ldsc_h2.tsv gencorr.sample_overlap.txt
    """
}
