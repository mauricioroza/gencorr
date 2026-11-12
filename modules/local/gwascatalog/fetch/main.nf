/*
 * GWASCATALOG_FETCH
 * -----------------
 * Downloads a study's harmonised summary statistics from the GWAS Catalog
 * (bin/fetch_gwascatalog.py): finds the file for the accession, downloads it
 * with its metadata, checks the md5 the Catalog publishes, and writes a JSON
 * with the provenance and the study's sample size, which HARMONIZE_SUMSTATS
 * records in the QC report (and uses as N when the samplesheet gives none).
 *
 * meta fields used: id, gwas_catalog (accession), ancestry (default EUR)
 * ext.args: e.g. '--base <url or folder>' (default: the GWAS Catalog FTP)
 * Needs network access. With --gwascatalog_cache, files are kept in
 * <cache>/<accession>/ and reused by later runs (storeDir).
 */
process GWASCATALOG_FETCH {
    tag "${meta.id}:${meta.gwas_catalog}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    val meta

    output:
    tuple val(meta), path("*.h.tsv.gz"), path("${meta.id}.gwascatalog.json"), emit: sumstats
    // no versions output: storeDir (--gwascatalog_cache) cannot be used with eval
    // outputs, and the script uses only Python's standard library (the same
    // Python as HARMONIZE_SUMSTATS, which reports its version)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fetch_gwascatalog.py \\
        --accession ${meta.gwas_catalog} \\
        --id ${meta.id} \\
        --ancestry ${meta.ancestry ?: 'EUR'} \\
        ${args}
    """

    stub:
    """
    printf 'rsid\\teffect_allele\\tother_allele\\tbeta\\tstandard_error\\tp_value\\n' | gzip > ${meta.gwas_catalog}.h.tsv.gz
    printf '{"source": "gwas_catalog", "accession": "${meta.gwas_catalog}", "n": 1000, "n_reason": null, "warnings": []}\\n' > ${meta.id}.gwascatalog.json
    """
}
