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
    SNP set: by default every harmonized SNP is used. With --lava_snp_list
    (e.g. w_hm3.snplist) the input is first restricted to that list, as
    BigaGWAS does with its HapMap3-based list.
    Results are identical however the loci are chunked: run_lava.R seeds each
    (locus, pair) separately (tests/python/test_lava_runner.py).
*/

include { EXPORT_SUMSTATS } from '../../../modules/local/export_sumstats'
include { LAVA_RUN        } from '../../../modules/local/lava/run'
include { MERGE_RESULTS   } from '../../../modules/local/merge_results'
include { traitAncestry   } from '../utils_nfcore_gencorr_pipeline'

workflow METHOD_LAVA {

    take:
    ch_harmonized    // channel: [ meta, harmonized ]
    ch_pairs         // channel: [ meta1, meta2 ]  (same-ancestry pairs only)
    ch_reference     // channel: [ ancestry, prefix, [ reference files ] ], one per ancestry
    ch_loci          // value channel: locus file (LOC CHR START STOP)
    ch_overlap       // value channel: sample-overlap matrix, or [] for none
    ch_snp_list      // value channel: SNP list to restrict LAVA input to, or [] for all SNPs
    snp_set          // string: 'pairwise' or 'joint'
    loci_per_chunk   // integer

    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'lava', ch_snp_list )

    // id -> LAVA-format file, as one map so groups can look up their files
    def ch_file_map = EXPORT_SUMSTATS.out.sumstats
        .map    { meta, f -> [ meta.id, f ] }
        .toList()
        .map    { list -> list.collectEntries() }

    // Trait groups: [ id, traits: [ metas ], pairs: [ [id1, id2], ... ], ancestry ]
    // Pairs are same-ancestry, so every group has one ancestry and one reference.
    def ch_groups
    if (snp_set == 'joint') {
        // one joint group per ancestry ('joint' when there is only one)
        ch_groups = ch_pairs
            .map        { m1, m2 -> [ traitAncestry(m1), [ m1, m2 ] ] }
            .groupTuple ( by: 0 )
            .toList     ()
            .flatMap    { per_ancestry ->
                per_ancestry.collect { anc, pairs ->
                    def traits = pairs.collectMany { pair -> pair }.unique { m -> m.id }.sort { m -> m.id }
                    [ id: per_ancestry.size() == 1 ? 'joint' : "joint_${anc}".toString(), traits: traits,
                      pairs: pairs.collect { pair -> [ pair[0].id, pair[1].id ] }, ancestry: anc ]
                }
            }
    } else {
        ch_groups = ch_pairs.map { m1, m2 ->
            [ id: "${m1.id}__${m2.id}", traits: [ m1, m2 ], pairs: [ [ m1.id, m2.id ] ], ancestry: traitAncestry(m1) ]
        }
    }

    // Locus chunks: loci.1.txt, loci.2.txt, ... each keeping the header line
    def ch_chunks = ch_loci
        .splitText( by: loci_per_chunk as int, keepHeader: true, file: true )

    def ch_tasks = ch_groups
        .combine( ch_file_map )
        .map    { group, fmap -> [ group, group.traits.collect { m -> fmap[m.id] } ] }
        .combine( ch_chunks )                                // [ group, files, chunk ]

    // each task with its ancestry's LD reference
    def ch_lava_in = ch_tasks
        .map      { group, files, chunk -> [ group.ancestry, group, files, chunk ] }
        .combine  ( ch_reference.map { anc, prefix, ref_files -> [ anc, [ prefix, ref_files ] ] }, by: 0 )
        .multiMap { _anc, group, files, chunk, ref ->
            task: [ group, files, chunk ]
            ref:  ref
        }

    LAVA_RUN ( ch_lava_in.task, ch_lava_in.ref, ch_overlap )

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
