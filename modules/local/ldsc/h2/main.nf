/*
 * LDSC_H2
 * -------
 * SNP-heritability of one trait. Always runs on the observed scale; for binary
 * traits with pop_prev and a sample prevalence (samp_prev, or n_cases and
 * n_controls) it runs a second time with
 * --samp-prev/--pop-prev, because LDSC then reports ONLY the liability-scale
 * value. Both numbers therefore come from LDSC itself.
 *
 * <id>.coverage.log (bin/ldsc_coverage.py): share of the LD score reference
 * the trait covers, used by the confidence tier.
 *
 * The h2 intercept from here is the diagonal of the sample-overlap matrix
 * passed to LAVA.
 *
 * ld: [ meta(prefix, split_by_chr), ref_ld_dir, w_ld_dir ]
 */
process LDSC_H2 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-ldsc:3.0.1-128f996"

    input:
    tuple val(meta), path(sumstats)
    tuple val(ld), path(ref_ld, stageAs: 'ref_ld'), path(w_ld, stageAs: 'w_ld')

    output:
    tuple val(meta), path("${meta.id}.h2*.log"), emit: log
    tuple val(meta), path("${meta.id}.coverage.log"), emit: coverage
    tuple val("${task.process}"), val('ldsc'), eval("grep -m1 -o 'Version [0-9.]*' ${meta.id}.h2.log | cut -d' ' -f2"), emit: versions_ldsc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def ref_opt = ld.split_by_chr ? '--ref-ld-chr' : '--ref-ld'
    def w_opt   = ld.split_by_chr ? '--w-ld-chr'   : '--w-ld'
    def ld_args = "${ref_opt} ref_ld/${ld.prefix} ${w_opt} w_ld/${ld.prefix}"
    // sample prevalence: samplesheet samp_prev (e.g. 0.5 for an effective N) wins
    // over n_cases / (n_cases + n_controls)
    def samp    = meta.samp_prev ? meta.samp_prev as double
                : (meta.n_cases && meta.n_controls ? (meta.n_cases as double) / ((meta.n_cases as double) + (meta.n_controls as double)) : null)
    def liab    = meta.pop_prev && samp != null
    """
    ldsc.py --h2 ${sumstats} ${ld_args} --out ${meta.id}.h2 ${args}

    ldsc_coverage.py --sumstats ${sumstats} --ref-ld ref_ld/${ld.prefix} ${ld.split_by_chr ? '--split-by-chr' : ''} --out ${meta.id}.coverage.log

    ${ liab ? "ldsc.py --h2 ${sumstats} ${ld_args} --samp-prev ${samp} --pop-prev ${meta.pop_prev} --out ${meta.id}.h2_liab ${args}" : '' }
    """

    stub:
    """
    printf '* Version 3.0.1\\n' > ${meta.id}.h2.log
    touch ${meta.id}.coverage.log
    """
}
