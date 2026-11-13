/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD_SUMHER - global genetic correlation with SumHer (LDAK)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Follows the gencorr method interface (docs/extending.md):

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]      (traits of the pairs)
             ch_pairs       [ meta1, meta2 ]                 (same-ancestry pairs)
             + method-specific reference inputs
      emit:  results        path: standard results table (one row per pair)

    Steps:  EXPORT_SUMSTATS(sumher) -> SUMHER_RG (one LDAK run per pair, with the
            tagging file of the pair's ancestry) -> SUMHER_SUMMARISE
*/

include { EXPORT_SUMSTATS  } from '../../../modules/local/export_sumstats'
include { SUMHER_RG        } from '../../../modules/local/sumher/rg'
include { SUMHER_SUMMARISE } from '../../../modules/local/sumher/summarise'
include { traitAncestry    } from '../utils_nfcore_gencorr_pipeline'

workflow METHOD_SUMHER {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]  (same-ancestry pairs only)
    ch_tagfile       // channel: [ ancestry, tagging file ], one per ancestry
    positions        // value channel: rsID -> GRCh37 position table (predictors named CHR:BP), or [] for rsIDs

    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'sumher', positions )

    def ch_by_id = EXPORT_SUMSTATS.out.sumstats.map { meta, f -> [ meta.id, f ] }

    // [ pair_meta, sumstats1, sumstats2 ] with the pair's ancestry
    def ch_pair_files = ch_pairs
        .map     { m1, m2 -> [ m1.id, m2.id, traitAncestry(m1) ] }
        .combine ( ch_by_id, by: 0 )                        // [ id1, id2, anc, f1 ]
        .map     { id1, id2, anc, f1 -> [ id2, id1, anc, f1 ] }
        .combine ( ch_by_id, by: 0 )                        // [ id2, id1, anc, f1, f2 ]
        .map     { id2, id1, anc, f1, f2 ->
            [ anc, [ id: "${id1}__${id2}".toString(), trait1: id1, trait2: id2, ancestry: anc ], f1, f2 ]
        }

    def ch_rg_in = ch_pair_files
        .combine  ( ch_tagfile, by: 0 )
        .multiMap { _anc, meta, f1, f2, tagfile ->
            pair:    [ meta, f1, f2 ]
            tagfile: tagfile
        }

    SUMHER_RG ( ch_rg_in.pair, ch_rg_in.tagfile )

    SUMHER_SUMMARISE (
        SUMHER_RG.out.log.map  { _meta, log  -> log  }.collect(),
        SUMHER_RG.out.cors.map { _meta, cors -> cors }.collect()
    )

    emit:
    results = SUMHER_SUMMARISE.out.results   // path: gencorr.sumher.tsv (standard results)
}
