/*
 * LDSC_RG
 * -------
 * Cross-trait LDSC for one trait1 against a batch of traits:
 *     ldsc.py --rg trait1,traitA,traitB,...
 * LDSC estimates each (trait1, traitX) pair independently inside one call, so
 * batching only saves start-up/IO cost; it does not change any estimate
 * (verified by tests/python/test_ldsc_batching.py).
 *
 * meta.id    : batch name, e.g. "bmi__b1"
 * ld         : [ meta(prefix, split_by_chr), ref_ld_dir, w_ld_dir ]
 */
process LDSC_RG {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(sumstats1), path(others)
    tuple val(ld), path(ref_ld, stageAs: 'ref_ld'), path(w_ld, stageAs: 'w_ld')

    output:
    tuple val(meta), path("${meta.id}.rg.log"), emit: log
    tuple val("${task.process}"), val('ldsc'), eval("grep -m1 -o 'Version [0-9.]*' ${meta.id}.rg.log | cut -d' ' -f2"), emit: versions_ldsc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def ref_opt = ld.split_by_chr ? '--ref-ld-chr' : '--ref-ld'
    def w_opt   = ld.split_by_chr ? '--w-ld-chr'   : '--w-ld'
    def files   = ([sumstats1] + (others instanceof List ? others : [others])).join(',')
    """
    ldsc.py \\
        --rg ${files} \\
        ${ref_opt} ref_ld/${ld.prefix} \\
        ${w_opt} w_ld/${ld.prefix} \\
        --out ${meta.id}.rg \\
        ${args}
    """

    stub:
    """
    printf '* Version 1.0.1\\n' > ${meta.id}.rg.log
    """
}
