/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD_POPCORN - cross-ancestry genetic correlation (Popcorn, Brown et al. 2016)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Follows the gencorr method interface (docs/extending.md), but takes the
    CROSS-ancestry pairs (the ones LDSC, LAVA and SumHer cannot estimate):

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]      (traits of the pairs)
             ch_pairs       [ meta1, meta2 ]                 (cross-ancestry pairs)
             + method-specific reference inputs
      emit:  results        path: standard results table (one row per pair)

    Each pair is ordered so trait 1 is from the score file's first population
    (cfile_ancestries, e.g. ['EUR', 'EAS'] for scores computed with
    --bfile1 EUR --bfile2 EAS).

    Steps:  EXPORT_SUMSTATS(popcorn) -> POPCORN_FIT (one run per pair)
            -> POPCORN_SUMMARISE
*/

include { EXPORT_SUMSTATS   } from '../../../modules/local/export_sumstats'
include { POPCORN_FIT       } from '../../../modules/local/popcorn/fit'
include { POPCORN_SUMMARISE } from '../../../modules/local/popcorn/summarise'
include { traitAncestry     } from '../utils_nfcore_gencorr_pipeline'

workflow METHOD_POPCORN {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]  (cross-ancestry pairs only)
    ch_cfile         // value channel: Popcorn score file for the two populations
    cfile_ancestries // list: ancestries of the score file's population 1 and 2, e.g. ['EUR', 'EAS']

    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'popcorn', [] )

    def ch_by_id = EXPORT_SUMSTATS.out.sumstats.map { meta, f -> [ meta.id, f ] }

    // pair meta, ordered to the score file; liability options for binary traits
    def ch_ordered = ch_pairs.map { m1, m2 ->
        def first  = traitAncestry(m1) == cfile_ancestries[0] ? m1 : m2
        def second = first.is(m1) ? m2 : m1
        def l1 = liabilityPrevalences(first)
        def l2 = liabilityPrevalences(second)
        def meta = [ id: "${first.id}__${second.id}".toString(), trait1: first.id, trait2: second.id,
                     ancestry_1: traitAncestry(first), ancestry_2: traitAncestry(second) ]
        if (l1) { meta.k1 = l1.k; meta.p1 = l1.p }
        if (l2) { meta.k2 = l2.k; meta.p2 = l2.p }
        [ first.id, second.id, meta ]
    }

    def ch_fit_in = ch_ordered
        .combine ( ch_by_id, by: 0 )                        // [ id1, id2, meta, f1 ]
        .map     { id1, id2, meta, f1 -> [ id2, id1, meta, f1 ] }
        .combine ( ch_by_id, by: 0 )                        // [ id2, id1, meta, f1, f2 ]
        .map     { _id2, _id1, meta, f1, f2 -> [ meta, f1, f2 ] }

    POPCORN_FIT ( ch_fit_in, ch_cfile )

    POPCORN_SUMMARISE (
        POPCORN_FIT.out.results.map { _meta, f    -> f    }.collect(),
        POPCORN_FIT.out.log.map     { _meta, logs -> logs }.flatten().collect()
    )

    emit:
    results = POPCORN_SUMMARISE.out.results   // path: gencorr.popcorn.tsv (standard results)
}

//
// Population (K) and sample (P) prevalence for Popcorn's liability scale, for a
// binary trait with n_cases, n_controls and pop_prev in the samplesheet; [:] otherwise.
//
def liabilityPrevalences(Map m) {
    if (!(m.n_cases && m.n_controls && m.pop_prev)) {
        return [:]
    }
    return [ k: m.pop_prev, p: (m.n_cases as double) / ((m.n_cases as double) + (m.n_controls as double)) ]
}
