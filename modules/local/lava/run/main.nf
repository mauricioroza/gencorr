/*
 * LAVA_RUN
 * --------
 * LAVA univariate + bivariate local rg for one chunk of loci and a set of
 * pairs (bin/run_lava.R). The subworkflow decides which traits share a task:
 *   lava_snp_set = 'pairwise' : one pair per task (SNPs shared by the two traits)
 *   lava_snp_set = 'joint'    : all traits per task (SNPs shared by all traits)
 *
 * meta        : [ id, traits: [ trait metas ], pairs: [[t1, t2], ...] ]
 * sumstats    : LAVA-format files, named <trait id>.lava.tsv
 * reference   : all files of the LD reference (PLINK bed/bim/fam or LAVA
 *               LD format), staged together; ref_prefix is their common prefix
 * loci        : locus file chunk (LOC CHR START STOP)
 * overlap     : sample-overlap matrix, or [] for none
 */
process LAVA_RUN {
    tag "${meta.id}:${loci.baseName}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/mauricioroza/gencorr-tools:0.1.0"

    input:
    tuple val(meta), path(sumstats, stageAs: 'sumstats/*'), path(loci)
    tuple val(ref_prefix), path(reference, stageAs: 'reference/*')
    path overlap

    output:
    tuple val(meta), path("${prefix}.lava.tsv")     , emit: results
    tuple val(meta), path("${prefix}.lava_univ.tsv"), emit: univ
    tuple val("${task.process}"), val('lava'), eval("Rscript -e 'cat(as.character(packageVersion(\"LAVA\")))'"), emit: versions_lava, topic: versions
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}.${loci.baseName}"
    // LAVA input.info: continuous traits get cases/controls = NA (LAVA treats
    // NA case proportion as continuous); 'prevalence' enables h2.latent.
    def info  = meta.traits.collect { t ->
        def cases = t.n_cases && t.n_controls ? t.n_cases    : 'NA'
        def ctrls = t.n_cases && t.n_controls ? t.n_controls : 'NA'
        def prev  = t.pop_prev ?: 'NA'
        "${t.id}\t${cases}\t${ctrls}\t${prev}\tsumstats/${t.id}.lava.tsv"
    }.join('\n')
    def pairs = meta.pairs.collect { p -> "${p[0]}\t${p[1]}" }.join('\n')
    def ov    = overlap ? "--overlap ${overlap}" : ''
    """
    printf 'phenotype\\tcases\\tcontrols\\tprevalence\\tfilename\\n${info}\\n' > input.info.txt
    printf 'trait1\\ttrait2\\n${pairs}\\n' > pairs.tsv

    run_lava.R \\
        --info input.info.txt \\
        --ref reference/${ref_prefix} \\
        --loci ${loci} \\
        --pairs pairs.tsv \\
        --out ${prefix} \\
        ${ov} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${loci.baseName}"
    """
    touch ${prefix}.lava.tsv ${prefix}.lava_univ.tsv
    """
}
