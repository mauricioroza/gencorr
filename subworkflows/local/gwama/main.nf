/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GWAMA - combine groups of correlated traits into new traits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    N-weighted multivariate GWAMA (Baselmans et al. 2019, docs/gwama.md).
    Not a METHOD_* subworkflow: it does not estimate rg, it makes new traits,
    which workflows/gencorr.nf adds to the harmonized traits before the methods
    run, so every method analyses them like samplesheet traits.

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]   (samplesheet traits)
             ch_groups      [ group_meta, [ member metas ] ]
             + LDSC references, as for METHOD_LDSC
      emit:  sumstats       [ group_meta, harmonized.tsv.gz ]  (combined traits)
             qc             [ group_meta, gwama.json ]

    Steps:  METHOD_LDSC on every pair inside a group (h2, univariate and
            cross-trait intercepts) -> GWAMA_COMBINE per group
*/

include { METHOD_LDSC as GWAMA_LDSC } from '../method_ldsc'
include { GWAMA_COMBINE             } from '../../../modules/local/gwama/combine'

workflow GWAMA {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_groups        // channel: [ group_meta, [ member metas ] ]
    ch_ld            // channel: [ ancestry, ld_meta, ref_ld_dir, w_ld_dir ], one per ancestry
    merge_alleles    // value channel: path to merge-alleles list, or [] for none
    rg_batch_size    // integer: max traits per ldsc --rg call

    main:
    // every pair of traits inside a group, once (groups may share traits)
    def ch_pairs = ch_groups
        .flatMap { _g, metas ->
            (0..<metas.size()).collectMany { i -> ((i + 1)..<metas.size()).collect { j -> [ metas[i], metas[j] ] } }
        }
        .unique  { pair -> [ pair[0].id, pair[1].id ].sort().join('\u0000') }

    // harmonized files of the traits that are in at least one group
    def ch_member_ids = ch_groups
        .flatMap { _g, metas -> metas.collect { m -> m.id } }
        .unique()
        .collect()
        .map     { ids -> [ ids ] }
    def ch_members = ch_harmonized
        .combine ( ch_member_ids )
        .filter  { meta, _f, ids -> meta.id in ids }
        .map     { meta, f, _ids -> [ meta, f ] }

    GWAMA_LDSC ( ch_members, ch_pairs, ch_ld, merge_alleles, rg_batch_size )

    // [ group_meta, [ files in the order of group_meta.gwama_traits ] ]
    def ch_files_by_id = ch_members.map { meta, f -> [ meta.id, f ] }
    def ch_group_files = ch_groups
        .flatMap { g, metas -> metas.withIndex().collect { m, i -> [ m.id, g, i ] } }
        .combine ( ch_files_by_id, by: 0 )                      // [ id, g, i, f ]
        .map     { _id, g, i, f -> [ groupKey(g, g.gwama_traits.size()), [ i, f ] ] }
        .groupTuple()
        .map     { g, items -> [ g.getGroupTarget(), items.sort { x -> x[0] }.collect { x -> x[1] } ] }

    GWAMA_COMBINE (
        ch_group_files,
        GWAMA_LDSC.out.h2.first(),
        GWAMA_LDSC.out.results.first()
    )

    emit:
    sumstats = GWAMA_COMBINE.out.sumstats   // [ group_meta, harmonized.tsv.gz ]
    qc       = GWAMA_COMBINE.out.qc         // [ group_meta, gwama.json ]
    ldsc_h2  = GWAMA_LDSC.out.h2            // path: LDSC h2 of the group traits
}
