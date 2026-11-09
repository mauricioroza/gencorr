/*
 * POPCORN_FIT
 * -----------
 * Cross-ancestry genetic correlation for one pair (Popcorn, Brown et al. 2016):
 *     popcorn fit -v 1 --cfile <scores> --sfile1 <pop1 trait> --sfile2 <pop2 trait> <pair>.popcorn
 * sumstats1 must be the trait of the score file's first population (the
 * subworkflow orders each pair accordingly). Options (--gen_effect, --use_mle)
 * come from ext.args (conf/modules.config). For binary traits with n_cases,
 * n_controls and pop_prev in the samplesheet, the liability-scale options
 * (--K1/--P1, --K2/--P2) are added from the pair meta.
 *
 * A failing Popcorn run does not stop the pipeline: an empty <pair>.popcorn is
 * left and POPCORN_SUMMARISE reports the pair as 'popcorn_error' with the error.
 *
 * meta    : [ id: '<trait1>__<trait2>', trait1, trait2, ancestry_1, ancestry_2,
 *             k1, p1, k2, p2 (liability; optional) ]
 * cfile   : Popcorn score file (plain or .gz) for the two populations
 */
process POPCORN_FIT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-popcorn:1.0.0-20a29e8"

    input:
    tuple val(meta), path(sumstats1), path(sumstats2)
    path cfile

    output:
    tuple val(meta), path("${meta.id}.popcorn")                , emit: results
    tuple val(meta), path("${meta.id}.popcorn.{o,e}.log")      , emit: log
    tuple val("${task.process}"), val('popcorn'), eval("grep -m1 -oE 'Popcorn version [0-9.]+' ${meta.id}.popcorn.o.log | cut -d' ' -f3"), emit: versions_popcorn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def gz    = cfile.name.endsWith('.gz')
    def cname = gz ? cfile.name[0..-4] : cfile.name
    def liab  = [ meta.k1 ? "--K1 ${meta.k1} --P1 ${meta.p1}" : '',
                  meta.k2 ? "--K2 ${meta.k2} --P2 ${meta.p2}" : '' ].join(' ')
    """
    ${gz ? "gunzip -c ${cfile} > ${cname}" : ''}
    if ! popcorn fit -v 1 \\
            --cfile ${cname} \\
            --sfile1 ${sumstats1} \\
            --sfile2 ${sumstats2} \\
            ${liab} ${args} \\
            ${meta.id}.popcorn > /dev/null 2>&1; then
        echo "Popcorn failed for ${meta.id}; reported as popcorn_error" >&2
        tail -n 3 ${meta.id}.popcorn.e.log >&2 || true
    fi
    touch ${meta.id}.popcorn ${meta.id}.popcorn.o.log ${meta.id}.popcorn.e.log
    """

    stub:
    """
    printf 'Popcorn version 1.1\\nInvoking command: popcorn fit --sfile1 ${sumstats1} --sfile2 ${sumstats2} ${meta.id}.popcorn\\n' > ${meta.id}.popcorn.o.log
    touch ${meta.id}.popcorn ${meta.id}.popcorn.e.log
    """
}
