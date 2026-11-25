/*
 * POPCORN_COMPUTE
 * ---------------
 * Popcorn score file (LD scores within and across two populations) from one
 * PLINK reference per population (Popcorn, Brown et al. 2016):
 *     popcorn compute -v 1 --bfile1 <pop1> --bfile2 <pop2> [--extract snps] <out>
 * Population 1 is meta.ancestry_1. Score files only depend on the two
 * references, the SNP list and the options, so they are computed once and kept
 * with --popcorn_cache (storeDir, conf/modules.config). The file name records
 * the ancestries and the model (impact/effect), so a cache never mixes models.
 * Options (--maf, --gen_effect; the latter must match POPCORN_FIT) come from
 * ext.args.
 *
 * meta    : [ id: '<anc1>_<anc2>', ancestry_1, ancestry_2, model: 'impact'|'effect' ]
 * prefix1 : PLINK prefix (basename) of the population-1 reference; files staged in pop1/
 * extract : SNP list (first whitespace-separated column used, e.g. w_hm3.snplist), or []
 */
process POPCORN_COMPUTE {
    tag "${meta.id}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8"

    input:
    tuple val(meta), val(prefix1), path(ref1, stageAs: 'pop1/*'), val(prefix2), path(ref2, stageAs: 'pop2/*')
    path extract

    output:
    tuple val(meta), path("${meta.id}.${meta.model}.cscore"), emit: cfile
    tuple val(meta), path("${meta.id}.${meta.model}.cscore.log"), emit: log
    tuple val("${task.process}"), val('popcorn'), eval("grep -m1 -oE 'Popcorn version [0-9.]+' ${meta.id}.${meta.model}.cscore.log | cut -d' ' -f3"), emit: versions_popcorn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out  = "${meta.id}.${meta.model}.cscore"
    """
    ${extract ? "zcat -f ${extract} | awk 'NF && \$1 != \"SNP\" { print \$1 }' > snps.extract" : ''}
    popcorn compute -v 1 \\
        --bfile1 pop1/${prefix1} \\
        --bfile2 pop2/${prefix2} \\
        ${extract ? '--extract snps.extract' : ''} \\
        ${args} \\
        ${out} > ${out}.log 2>&1 || { tail -n 20 ${out}.log >&2; exit 1; }
    """

    stub:
    def out = "${meta.id}.${meta.model}.cscore"
    """
    printf 'Popcorn version 1.1\\n' > ${out}.log
    touch ${out}
    """
}
