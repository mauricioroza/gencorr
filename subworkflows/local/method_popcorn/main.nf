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

    Every ancestry pair needs a score file. Given ones arrive in ch_cfiles; the
    others are computed here (POPCORN_COMPUTE) from one PLINK reference per
    population. Each trait pair is ordered so trait 1 is from its score file's
    first population.

    Steps:  [POPCORN_COMPUTE] ; EXPORT_SUMSTATS(popcorn) -> POPCORN_FIT (one run
            per pair, with its ancestry pair's score file) -> POPCORN_SUMMARISE
*/

include { EXPORT_SUMSTATS   } from '../../../modules/local/export_sumstats'
include { POPCORN_COMPUTE   } from '../../../modules/local/popcorn/compute'
include { POPCORN_FIT       } from '../../../modules/local/popcorn/fit'
include { POPCORN_SUMMARISE } from '../../../modules/local/popcorn/summarise'
include { traitAncestry     } from '../utils_nfcore_gencorr_pipeline'

workflow METHOD_POPCORN {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]  (cross-ancestry pairs only)
    ch_cfiles        // channel: [ [ ancestry_1, ancestry_2 ], score file ]  (given score files, file order)
    ch_compute       // channel: [ meta(id, ancestry_1, ancestry_2, model), prefix1, [ ref1 ], prefix2, [ ref2 ] ]
    ch_extract       // value channel: SNP list for computed scores, or []

    main:
    POPCORN_COMPUTE ( ch_compute, ch_extract )

    // [ key, [ pop1, pop2 ], score file ]; key = the two ancestries in any order
    def ch_cfile_by_key = ch_cfiles
        .mix( POPCORN_COMPUTE.out.cfile.map { meta, f -> [ [ meta.ancestry_1, meta.ancestry_2 ], f ] } )
        .map { pops, f -> [ ancestryKey(pops[0], pops[1]), pops, f ] }

    EXPORT_SUMSTATS ( ch_harmonized, 'popcorn', [] )

    def ch_by_id = EXPORT_SUMSTATS.out.sumstats.map { meta, f -> [ meta.id, f ] }

    // pair meta, ordered to its score file; liability options for binary traits
    def ch_ordered = ch_pairs
        .map     { m1, m2 -> [ ancestryKey(traitAncestry(m1), traitAncestry(m2)), m1, m2 ] }
        .combine ( ch_cfile_by_key, by: 0 )                 // [ key, m1, m2, pops, cfile ]
        .map     { _key, m1, m2, pops, cfile ->
            def first  = traitAncestry(m1) == pops[0] ? m1 : m2
            def second = first.is(m1) ? m2 : m1
            def l1 = liabilityPrevalences(first)
            def l2 = liabilityPrevalences(second)
            def meta = [ id: "${first.id}__${second.id}".toString(), trait1: first.id, trait2: second.id,
                         ancestry_1: traitAncestry(first), ancestry_2: traitAncestry(second) ]
            if (l1) { meta.k1 = l1.k; meta.p1 = l1.p }
            if (l2) { meta.k2 = l2.k; meta.p2 = l2.p }
            [ first.id, second.id, meta, cfile ]
        }

    def ch_fit_in = ch_ordered
        .combine ( ch_by_id, by: 0 )                        // [ id1, id2, meta, cfile, f1 ]
        .map     { id1, id2, meta, cfile, f1 -> [ id2, id1, meta, cfile, f1 ] }
        .combine ( ch_by_id, by: 0 )                        // [ id2, id1, meta, cfile, f1, f2 ]
        .map     { _id2, _id1, meta, cfile, f1, f2 -> [ meta, f1, f2, cfile ] }

    POPCORN_FIT ( ch_fit_in )

    POPCORN_SUMMARISE (
        POPCORN_FIT.out.results.map { _meta, f    -> f    }.collect(),
        POPCORN_FIT.out.log.map     { _meta, logs -> logs }.flatten().collect()
    )

    emit:
    results = POPCORN_SUMMARISE.out.results   // path: gencorr.popcorn.tsv (standard results)
}

// Order-independent key of an ancestry pair, e.g. 'EAS|EUR'
def ancestryKey(String a, String b) {
    return [ a, b ].sort().join('|')
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
