/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD_LAVA - local genetic correlation (LAVA, Werme et al. 2022)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Follows the gencorr method interface (docs/extending.md):

      take:  ch_harmonized  [ meta, harmonized.tsv.gz ]
             ch_pairs       [ meta1, meta2 ]
             + method-specific reference inputs
      emit:  results        path: standard results table (one row per pair x locus)

    Parallelism: the locus file is split into chunks of `loci_per_chunk` loci,
    and tasks are formed as (trait group x chunk):
      snp_set = 'pairwise' : one group per pair. Each pair uses the SNPs shared by
                             its two traits, so its result does not depend on which
                             other traits are in the run. (default)
      snp_set = 'joint'    : one group with every trait. All pairs use the SNPs
                             shared by ALL traits (LAVA's behaviour when many
                             phenotypes are loaded at once, as in BigaGWAS batches).
    Results are identical however the loci are chunked: run_lava.R seeds each
    (locus, pair) separately (tests/python/test_lava_runner.py).
*/

include { EXPORT_SUMSTATS } from '../../../modules/local/export_sumstats'
include { LAVA_RUN        } from '../../../modules/local/lava/run'
include { MERGE_RESULTS   } from '../../../modules/local/merge_results'

workflow METHOD_LAVA {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]
    ch_reference     // value channel: [ prefix, [ reference files ] ]
    ch_loci          // value channel: locus file (LOC CHR START STOP)
    ch_overlap       // value channel: sample-overlap matrix, or [] for none
    snp_set          // string: 'pairwise' or 'joint'
    loci_per_chunk   // integer

    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'lava' )

    // id -> LAVA-format file, as one map so groups can look up their files
    def ch_file_map = EXPORT_SUMSTATS.out.sumstats
        .map    { meta, f -> [ meta.id, f ] }
        .toList()
        .map    { list -> list.collectEntries() }

    // Trait groups: [ id, traits: [ metas ], pairs: [ [id1, id2], ... ] ]
    def ch_groups
    if (snp_set == 'joint') {
        ch_groups = ch_pairs
            .toList()
            .map { pairs ->
                def traits = pairs.collectMany { pair -> pair }.unique { m -> m.id }.sort { m -> m.id }
                [ id: 'joint', traits: traits, pairs: pairs.collect { pair -> [ pair[0].id, pair[1].id ] } ]
            }
    } else {
        ch_groups = ch_pairs.map { m1, m2 ->
            [ id: "${m1.id}__${m2.id}", traits: [ m1, m2 ], pairs: [ [ m1.id, m2.id ] ] ]
        }
    }

    // Locus chunks: loci.1.txt, loci.2.txt, ... each keeping the header line
    def ch_chunks = ch_loci
        .splitText( by: loci_per_chunk as int, keepHeader: true, file: true )

    def ch_tasks = ch_groups
        .combine( ch_file_map )
        .map    { group, fmap -> [ group, group.traits.collect { m -> fmap[m.id] } ] }
        .combine( ch_chunks )                                // [ group, files, chunk ]

    LAVA_RUN ( ch_tasks, ch_reference, ch_overlap )

    MERGE_RESULTS (
        LAVA_RUN.out.results
            .map { _meta, table -> table }
            .collect( sort: { a, b -> a.name <=> b.name } )
            .map { tables -> [ [ id: 'gencorr.lava' ], tables ] }
    )

    def ch_univ = LAVA_RUN.out.univ
        .map { _meta, table -> table }
        .collectFile( name: 'gencorr.lava_univ.tsv', keepHeader: true, skip: 1, sort: { f -> f.name },
                      storeDir: "${params.outdir}/lava" )

    emit:
    results = MERGE_RESULTS.out.results.map { _meta, f -> f }   // path: standard results
    univ    = ch_univ                                          // path: local univariate h2 per locus
}
