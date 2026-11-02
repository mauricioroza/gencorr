//
// Subworkflow with functionality specific to the mauricioroza/gencorr pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD REGISTRY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Every method the pipeline can run, with the parameters it needs.
    To add a method: add an entry here, a subworkflow under subworkflows/local/,
    and one line in workflows/gencorr.nf. See docs/extending.md.
*/
def methodRegistry() {
    return [
        ldsc: [ scope: 'global', required_params: ['ldsc_ld'],             conda_ok: true  ],
        lava: [ scope: 'local',  required_params: ['lava_ref', 'lava_loci'], conda_ok: false ],
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def before_text = ""
    def after_text = ""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR> --methods ldsc,lava"
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Pipeline-specific parameter checks (things JSON schema cannot express)
    //
    def methods = validateGencorrParams(params, workflow.profile)   // the workflow re-reads params.methods

    //
    // Samplesheet -> [ meta, sumstats ]. Empty optional columns become absent
    // meta keys, so modules can test them with `meta.n_cases ? ... : ...`.
    //
    def samples = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, sumstats -> [ cleanMeta(meta), sumstats ] }
    validateSamples(samples)

    //
    // Trait pairs, decided once here so every method analyses the same pairs.
    //
    def ids   = samples.collect { it[0].id }
    def pair_list = makePairs(ids, params.pair_mode, params.targets,
                              params.pairs ? readPairsFile(params.pairs) : null)
    log.info "gencorr: ${samples.size()} traits, ${pair_list.size()} pairs (pair_mode = ${params.pair_mode}), methods: ${methods.join(', ')}"

    def meta_by_id = samples.collectEntries { [ (it[0].id): it[0] ] }

    emit:
    samplesheet = channel.fromList(samples)                                        // [ meta, sumstats ]
    pairs       = channel.fromList(pair_list.collect { p -> [ meta_by_id[p[0]], meta_by_id[p[1]] ] }) // [ meta1, meta2 ]
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    workflow.onComplete {
        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. See the .command.err / .command.log of the failed task (path printed above)."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pure functions (no channels), unit-tested in
    subworkflows/local/utils_nfcore_gencorr_pipeline/tests/main.function.nf.test
*/

//
// Drop null / empty-list meta values produced by empty samplesheet cells
//
def cleanMeta(Map meta) {
    return meta.findAll { k, v -> v != null && v != [] && v != '' }
}

//
// Samplesheet consistency checks that the JSON schema cannot express
//
def validateSamples(List samples) {
    samples.each { meta, sumstats ->
        if (meta.n_cases && meta.n_controls && meta.n && Math.abs((meta.n as double) - (meta.n_cases as double) - (meta.n_controls as double)) > 0.5) {
            error("Samplesheet row '${meta.id}': n (${meta.n}) differs from n_cases + n_controls (${(meta.n_cases as double) + (meta.n_controls as double)}). Leave n empty for binary traits.")
        }
    }
    if (samples.size() < 2) {
        error("The samplesheet must contain at least 2 traits (found ${samples.size()}).")
    }
}

//
// Build the list of trait pairs.
//   all    : every unordered pair, in samplesheet order (i < j)
//   target : each target trait against every other trait (BigaGWAS "one vs group")
//   pairs  : exactly the pairs listed in a CSV (columns trait1,trait2)
// Pairs are unordered: [a, b] and [b, a] are the same pair; the first seen wins.
//
def makePairs(List ids, String mode, String targets, List listed) {
    def out  = []
    def seen = [] as Set
    def add  = { a, b ->
        if (a == b) error("Pair '${a}' vs '${b}': a trait cannot be paired with itself.")
        def key = [a, b].sort().join('\u0000')
        if (!seen.contains(key)) { seen << key; out << [a, b] }
    }
    switch (mode) {
        case 'all':
            for (int i = 0; i < ids.size(); i++)
                for (int j = i + 1; j < ids.size(); j++)
                    add(ids[i], ids[j])
            break
        case 'target':
            if (!targets) error("--pair_mode target requires --targets (comma-separated trait IDs).")
            def tlist = targets.tokenize(',')*.trim()
            def unknown = tlist.findAll { !(it in ids) }
            if (unknown) error("--targets: unknown trait ID(s) ${unknown}. Known: ${ids}")
            tlist.each { t -> ids.each { o -> if (o != t) add(t, o) } }
            break
        case 'pairs':
            if (listed == null) error("--pair_mode pairs requires --pairs <csv with trait1,trait2>.")
            listed.each { p ->
                def unknown = p.findAll { !(it in ids) }
                if (unknown) error("--pairs: unknown trait ID(s) ${unknown} in pair ${p}. Known: ${ids}")
                add(p[0], p[1])
            }
            break
        default:
            error("--pair_mode must be one of all, target, pairs (got '${mode}').")
    }
    if (!out) error("No trait pairs to analyse with pair_mode '${mode}'.")
    return out
}

def readPairsFile(String path) {
    def lines = file(path, checkIfExists: true).readLines().findAll { it.trim() }
    def header = lines[0].split(/[,\t]/)*.trim()
    if (header[0..1] != ['trait1', 'trait2']) error("--pairs: header must start with trait1,trait2 (got ${header})")
    return lines.drop(1).collect { it.split(/[,\t]/)*.trim()[0..1] }
}

//
// Cross-parameter checks. Returns the list of methods to run.
//
def validateGencorrParams(Map p, String profile) {
    def registry = methodRegistry()
    def methods  = (p.methods ?: '').tokenize(',')*.trim().unique()
    def unknown  = methods.findAll { !registry.containsKey(it) }
    if (!methods) error("--methods is empty. Choose from: ${registry.keySet().join(', ')}")
    if (unknown)  error("--methods: unknown method(s) ${unknown}. Choose from: ${registry.keySet().join(', ')}")

    methods.each { m ->
        registry[m].required_params.each { rp ->
            if (!p[rp]) error("--methods ${m} requires --${rp}.")
        }
        def conda = (profile ?: '').tokenize(',').intersect(['conda', 'mamba'])
        if (conda && !registry[m].conda_ok) {
            error("Method '${m}' has no conda package (LAVA is installed from GitHub). Use -profile docker, singularity or apptainer.")
        }
    }
    if ('lava' in methods && p.lava_overlap == 'ldsc' && !('ldsc' in methods)) {
        error("--lava_overlap ldsc estimates sample overlap from LDSC intercepts, so it needs --methods to include ldsc. Alternatively use --lava_overlap none or a matrix file.")
    }
    if (!(p.lava_snp_set in ['pairwise', 'joint'])) {
        error("--lava_snp_set must be 'pairwise' or 'joint' (got '${p.lava_snp_set}').")
    }
    if (p.ldsc_ld && !p.ldsc_merge_alleles && 'ldsc' in methods) {
        log.warn "No --ldsc_merge_alleles given. LDSC fails a pair when an rsID has different alleles in the two traits; the usual choice is w_hm3.snplist."
    }
    return methods
}

//
// Citations for the tools actually used (printed at the end of docs/output.md)
//
def toolCitationText() {
    return [
        "Tools used in the workflow included:",
        "LDSC (Bulik-Sullivan et al. 2015),",
        "LAVA (Werme et al. 2022)."
    ].join(' ').trim()
}
