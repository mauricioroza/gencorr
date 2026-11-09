/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENCORR - main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. HARMONIZE_SUMSTATS  every trait -> one harmonized format (+ QC report)
    2. METHOD_<X>          each selected method, all on the same harmonized
                           files and the same trait pairs
    3. MERGE_RESULTS       all methods -> results/gencorr_results.tsv
                           (one standard table, see docs/output.md)

    Adding a method = a subworkflow with the METHOD_* interface + one block
    below + an entry in methodRegistry() (docs/extending.md).
*/

include { HARMONIZE_SUMSTATS          } from '../modules/local/harmonize_sumstats'
include { MERGE_RESULTS as MERGE_ALL  } from '../modules/local/merge_results'
include { METHOD_LDSC                 } from '../subworkflows/local/method_ldsc'
include { METHOD_LAVA                 } from '../subworkflows/local/method_lava'
include { METHOD_SUMHER               } from '../subworkflows/local/method_sumher'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { traitAncestry; ancestryParam; methodRegistry } from '../subworkflows/local/utils_nfcore_gencorr_pipeline'

workflow GENCORR {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    ch_pairs         // channel: [ meta1, meta2 ]
    outdir

    main:
    def methods     = params.methods.tokenize(',').collect { m -> m.trim() }
    def ch_versions = channel.empty()
    def ch_results  = channel.empty()
    def ch_overlap  = channel.value([])     // LAVA sample overlap: none unless set below

    //
    // 1. Harmonize all summary statistics
    //
    //    Traits whose IDs are not rsIDs get them from the rsID map for their
    //    build (samplesheet column 'build', default GRCh37). No liftover: every
    //    method matches variants by rsID and uses its reference's positions.
    def rsid_maps = [
        GRCh37: params.rsid_map_grch37 ? file(params.rsid_map_grch37, checkIfExists: true) : [],
        GRCh38: params.rsid_map_grch38 ? file(params.rsid_map_grch38, checkIfExists: true) : [],
    ]
    def ch_harmonize_in = ch_samplesheet.map { meta, sumstats ->
        [ meta, sumstats, rsid_maps[meta.build ?: 'GRCh37'] ]
    }
    HARMONIZE_SUMSTATS ( ch_harmonize_in )
    def ch_harmonized = HARMONIZE_SUMSTATS.out.sumstats

    //
    // Ancestry (samplesheet column, default EUR): LDSC and LAVA assume both traits
    // come from one population, so they get same-ancestry pairs, each with its
    // ancestry's references. Cross-ancestry pairs are kept in the results with
    // status 'cross_ancestry' (a cross-ancestry method such as Popcorn is needed).
    //
    def ch_pairs_split = ch_pairs.branch { m1, m2 ->
        same:  traitAncestry(m1) == traitAncestry(m2)
        cross: true
    }
    def ch_same_pairs = ch_pairs_split.same
    ch_results = ch_results.mix(
        ch_pairs_split.cross
            .flatMap { m1, m2 -> methods.collect { m -> crossAncestryRow(m, m1, m2) } }
            .collectFile( name: 'gencorr.cross_ancestry.tsv', seed: resultHeader() + '\n',
                          newLine: true, sort: true )
    )

    //
    // 2a. LDSC (global)
    //
    if ('ldsc' in methods) {
        // one entry per ancestry with LD scores: [ ancestry, ld_meta, ref_ld_dir, w_ld_dir ]
        def ch_ld = channel.fromList( ['EUR', 'EAS'].findAll { a -> params[ancestryParam('ldsc_ld', a)] }.collect { a ->
            [ a,
              [ prefix: params.ldsc_ld_prefix ?: '', split_by_chr: params.ldsc_ld_split_by_chr ],
              file(params[ancestryParam('ldsc_ld', a)], checkIfExists: true),
              file(params[ancestryParam('ldsc_w_ld', a)] ?: params[ancestryParam('ldsc_ld', a)], checkIfExists: true) ]
        } )
        def ch_merge_alleles = params.ldsc_merge_alleles
            ? channel.value(file(params.ldsc_merge_alleles, checkIfExists: true))
            : channel.value([])

        METHOD_LDSC ( ch_harmonized, ch_same_pairs, ch_ld, ch_merge_alleles, params.ldsc_rg_batch_size )
        ch_results = ch_results.mix( METHOD_LDSC.out.results )

        if (params.lava_overlap == 'ldsc') {
            ch_overlap = METHOD_LDSC.out.overlap
        }
    }
    if (params.lava_overlap && !(params.lava_overlap in ['ldsc', 'none'])) {
        ch_overlap = channel.value( file(params.lava_overlap, checkIfExists: true) )
    }

    //
    // 2b. LAVA (local)
    //
    if ('lava' in methods) {
        // one entry per ancestry with a reference: [ ancestry, prefix, [ files ] ]
        def ch_reference = channel.fromList( ['EUR', 'EAS'].findAll { a -> params[ancestryParam('lava_ref', a)] }.collect { a ->
            [ a ] + resolveLavaReference(params[ancestryParam('lava_ref', a)])
        } )
        def ch_loci      = channel.value( file(params.lava_loci, checkIfExists: true) )
        def ch_snp_list  = params.lava_snp_list
            ? channel.value( file(params.lava_snp_list, checkIfExists: true) )
            : channel.value([])

        METHOD_LAVA (
            ch_harmonized, ch_same_pairs, ch_reference, ch_loci, ch_overlap, ch_snp_list,
            params.lava_snp_set, params.lava_loci_per_chunk
        )
        ch_results = ch_results.mix( METHOD_LAVA.out.results )
    }

    //
    // 2c. SumHer (global, LDAK)
    //
    if ('sumher' in methods) {
        // one entry per ancestry with a tagging file: [ ancestry, tagfile ]
        def ch_tagfile = channel.fromList( ['EUR', 'EAS'].findAll { a -> params[ancestryParam('sumher_tagfile', a)] }.collect { a ->
            [ a, file(params[ancestryParam('sumher_tagfile', a)], checkIfExists: true) ]
        } )
        def ch_positions = params.sumher_snp_positions
            ? channel.value( file(params.sumher_snp_positions, checkIfExists: true) )
            : channel.value([])

        METHOD_SUMHER ( ch_harmonized, ch_same_pairs, ch_tagfile, ch_positions )
        ch_results = ch_results.mix( METHOD_SUMHER.out.results )
    }

    //
    // 3. One table for all methods
    //
    MERGE_ALL (
        ch_results
            .collect( sort: { a, b -> a.name <=> b.name } )
            .map { tables -> [ [ id: 'gencorr_results' ], tables ] }
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'gencorr_software_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    results  = MERGE_ALL.out.results   // [ meta, gencorr_results.tsv ]
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Header of the standard results table (bin/gencorr_common.py RESULT_COLUMNS;
// MERGE_RESULTS validates every table against it, so a mismatch fails loudly).
//
def resultHeader() {
    return ['method', 'scope', 'trait1', 'trait2', 'locus', 'chr', 'start', 'stop', 'n_snps',
            'rg', 'rg_se', 'rg_ci_low', 'rg_ci_high', 'z', 'p', 'h2_1', 'h2_1_se', 'h2_2', 'h2_2_se',
            'status', 'extra'].join('\t')
}

//
// One results row for a pair the method cannot estimate because the two traits
// have different ancestry.
//
def crossAncestryRow(String method, Map m1, Map m2) {
    def extra = "{\"ancestry_1\": \"${traitAncestry(m1)}\", \"ancestry_2\": \"${traitAncestry(m2)}\"}"
    def cells = [ method, methodRegistry()[method].scope, m1.id, m2.id ] + (['NA'] * 15) + [ 'cross_ancestry', extra ]
    return cells.join('\t')
}

//
// LAVA reference = a prefix, e.g. /refs/g1000_eur (PLINK) or /refs/lava-ukb-v1.1
// (LAVA LD format, per-chromosome files). All files starting with the prefix are
// staged together. For http(s) prefixes only PLINK .bed/.bim/.fam are fetched.
// Returns [ prefix_basename, [ files ] ].
//
def resolveLavaReference(String prefix) {
    def ref_files
    if (prefix ==~ /^(https?|ftp):\/\/.*/) {
        ref_files = ['bed', 'bim', 'fam'].collect { ext -> file("${prefix}.${ext}") }
    } else {
        ref_files = files("${prefix}*").findAll { f -> f.isFile() }
    }
    if (!ref_files) {
        error("--lava_ref: no reference files found with prefix '${prefix}' (expected <prefix>.bed/.bim/.fam or LAVA LD-format files).")
    }
    return [ file(prefix).name, ref_files ]
}
