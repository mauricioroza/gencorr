/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD_LDSC - global genetic correlation with LD score regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Follows the gencorr method interface (docs/extending.md):

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]      (all traits)
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

workflow METHOD_LDSC {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]
    ch_ld            // value channel: [ ld_meta(prefix, split_by_chr), ref_ld_dir, w_ld_dir ]
    merge_alleles    // value channel: path to merge-alleles list, or [] for none
    rg_batch_size    // integer: max traits per ldsc --rg call

    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'ldsc', [] )   // LDSC restricts SNPs via --merge-alleles in munge
    LDSC_MUNGE      ( EXPORT_SUMSTATS.out.sumstats, merge_alleles )
    LDSC_H2         ( LDSC_MUNGE.out.sumstats, ch_ld )

    //
    // Batch pairs by trait1: [ batch_meta, trait1.sumstats.gz, [ others ] ]
    // Batch names and contents are deterministic (others sorted by ID), so
    // -resume and reruns reuse the same tasks.
    //
    def ch_munged_by_id = LDSC_MUNGE.out.sumstats.map { meta, f -> [ meta.id, f ] }

    def ch_batches = ch_pairs
        .map     { m1, m2 -> [ m2.id, m1.id ] }
        .combine ( ch_munged_by_id, by: 0 )                 // [ id2, id1, f2 ]
        .map     { id2, id1, f2 -> [ id1, id2, f2 ] }
        .combine ( ch_munged_by_id, by: 0 )                 // [ id1, id2, f2, f1 ]
        .map     { id1, id2, f2, f1 -> [ id1, f1, id2, f2 ] }
        .groupTuple( by: 0 )                                // [ id1, [f1..], [id2..], [f2..] ]
        .flatMap { id1, f1s, id2s, f2s ->
            def others = [ id2s, f2s ].transpose().sort { x -> x[0] }
            others.collate( rg_batch_size as int ).withIndex().collect { batch, i ->
                [ [ id: "${id1}__b${i + 1}", trait1: id1, traits2: batch.collect { x -> x[0] } ], f1s[0], batch.collect { x -> x[1] } ]
            }
        }

    LDSC_RG ( ch_batches, ch_ld )

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
