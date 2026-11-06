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
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'

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
    HARMONIZE_SUMSTATS ( ch_samplesheet )
    def ch_harmonized = HARMONIZE_SUMSTATS.out.sumstats

    //
    // 2a. LDSC (global)
    //
    if ('ldsc' in methods) {
        def ch_ld = channel.value([
            [ prefix: params.ldsc_ld_prefix ?: '', split_by_chr: params.ldsc_ld_split_by_chr ],
            file(params.ldsc_ld, checkIfExists: true),
            file(params.ldsc_w_ld ?: params.ldsc_ld, checkIfExists: true)
        ])
        def ch_merge_alleles = params.ldsc_merge_alleles
            ? channel.value(file(params.ldsc_merge_alleles, checkIfExists: true))
            : channel.value([])

        METHOD_LDSC ( ch_harmonized, ch_pairs, ch_ld, ch_merge_alleles, params.ldsc_rg_batch_size )
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
        def ch_reference = channel.value( resolveLavaReference(params.lava_ref) )
        def ch_loci      = channel.value( file(params.lava_loci, checkIfExists: true) )
        def ch_snp_list  = params.lava_snp_list
            ? channel.value( file(params.lava_snp_list, checkIfExists: true) )
            : channel.value([])

        METHOD_LAVA (
            ch_harmonized, ch_pairs, ch_reference, ch_loci, ch_overlap, ch_snp_list,
            params.lava_snp_set, params.lava_loci_per_chunk
        )
        ch_results = ch_results.mix( METHOD_LAVA.out.results )
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
