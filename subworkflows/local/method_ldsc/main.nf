/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD_LDSC - global genetic correlation with LD score regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Follows the gencorr method interface (docs/extending.md):

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]      (traits of the pairs)
             ch_pairs       [ meta1, meta2 ]                 (pairs to analyse)
             + method-specific reference inputs
      emit:  results        path: standard results table (one row per pair)

    Steps:  EXPORT_SUMSTATS(ldsc) -> LDSC_MUNGE -> LDSC_H2 (per trait)
                                                -> LDSC_RG (batched per trait1)
            -> LDSC_SUMMARISE (results, per-trait h2, sample-overlap matrix)
*/

include { EXPORT_SUMSTATS } from '../../../modules/local/export_sumstats'
include { LDSC_MUNGE      } from '../../../modules/local/ldsc/munge'
include { LDSC_H2         } from '../../../modules/local/ldsc/h2'
include { LDSC_RG         } from '../../../modules/local/ldsc/rg'
include { LDSC_SUMMARISE  } from '../../../modules/local/ldsc/summarise'
include { traitAncestry   } from '../utils_nfcore_gencorr_pipeline'

workflow METHOD_LDSC {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]  (same-ancestry pairs only)
    ch_ld            // channel: [ ancestry, ld_meta(prefix, split_by_chr), ref_ld_dir, w_ld_dir ], one per ancestry
    merge_alleles    // value channel: path to merge-alleles list, or [] for none
    rg_batch_size    // integer: max traits per ldsc --rg call

    main:
    def ch_ld_by_ancestry = ch_ld.map { anc, ld_meta, ref_ld, w_ld -> [ anc, [ ld_meta, ref_ld, w_ld ] ] }

    EXPORT_SUMSTATS ( ch_harmonized, 'ldsc', [] )   // LDSC restricts SNPs via --merge-alleles in munge
    LDSC_MUNGE      ( EXPORT_SUMSTATS.out.sumstats, merge_alleles )

    // h2 per trait with its ancestry's LD scores (ch_harmonized holds only traits
    // of same-ancestry pairs, whose LD scores are checked at start-up)
    def ch_h2_in = LDSC_MUNGE.out.sumstats
        .map      { meta, f -> [ traitAncestry(meta), meta, f ] }
        .combine  ( ch_ld_by_ancestry, by: 0 )
        .multiMap { _anc, meta, f, ld ->
            sumstats: [ meta, f ]
            ld:       ld
        }
    LDSC_H2 ( ch_h2_in.sumstats, ch_h2_in.ld )

    //
    // Batch pairs by trait1: [ batch_meta, trait1.sumstats.gz, [ others ] ]
    // Batch names and contents are deterministic (others sorted by ID), so
    // -resume and reruns reuse the same tasks.
    //
    def ch_munged_by_id = LDSC_MUNGE.out.sumstats.map { meta, f -> [ meta.id, f ] }

    def ch_batches = ch_pairs
        .map     { m1, m2 -> [ m2.id, m1.id, traitAncestry(m1) ] }
        .combine ( ch_munged_by_id, by: 0 )                 // [ id2, id1, anc, f2 ]
        .map     { id2, id1, anc, f2 -> [ id1, id2, anc, f2 ] }
        .combine ( ch_munged_by_id, by: 0 )                 // [ id1, id2, anc, f2, f1 ]
        .map     { id1, id2, anc, f2, f1 -> [ id1, f1, id2, f2, anc ] }
        .groupTuple( by: 0 )                                // [ id1, [f1..], [id2..], [f2..], [anc..] ]
        .flatMap { id1, f1s, id2s, f2s, ancs ->
            def others = [ id2s, f2s ].transpose().sort { x -> x[0] }
            others.collate( rg_batch_size as int ).withIndex().collect { batch, i ->
                [ [ id: "${id1}__b${i + 1}", trait1: id1, traits2: batch.collect { x -> x[0] }, ancestry: ancs[0] ],
                  f1s[0], batch.collect { x -> x[1] } ]
            }
        }

    // each batch with its ancestry's LD scores (pairs are same-ancestry)
    def ch_rg_in = ch_batches
        .map      { bmeta, f1, others -> [ bmeta.ancestry, bmeta, f1, others ] }
        .combine  ( ch_ld_by_ancestry, by: 0 )
        .multiMap { _anc, bmeta, f1, others, ld ->
            batch: [ bmeta, f1, others ]
            ld:    ld
        }

    LDSC_RG ( ch_rg_in.batch, ch_rg_in.ld )

    LDSC_SUMMARISE (
        LDSC_H2.out.log.map { _meta, logs -> logs }.flatten().collect(),
        LDSC_RG.out.log.map { _meta, log -> log }.collect().ifEmpty([])
    )

    emit:
    results    = LDSC_SUMMARISE.out.results   // path: gencorr.ldsc.tsv (standard results)
    h2         = LDSC_SUMMARISE.out.h2        // path: per-trait LDSC h2
    overlap    = LDSC_SUMMARISE.out.overlap   // path: sample-overlap matrix for LAVA
    munge_logs = LDSC_MUNGE.out.log           // [ meta, log ]
}
