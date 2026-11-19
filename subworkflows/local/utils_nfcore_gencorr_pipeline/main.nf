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
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD REGISTRY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Every method the pipeline can run, with the parameters it needs.
    ancestry_params are per-ancestry references: the plain name is the EUR
    reference, '<name>_eas' the EAS one (see ancestryParam()). They are required
    only for ancestries that have at least one same-ancestry pair.
    conda_ok: false stops the run at start-up under -profile conda/mamba.
    cross_ancestry: true = the method takes the cross-ancestry pairs instead of
    the same-ancestry ones.
    To add a method: add an entry here, a subworkflow under subworkflows/local/,
    and one block in workflows/gencorr.nf. See docs/extending.md.
*/
def methodRegistry() {
    return [
        ldsc: [ scope: 'global', required_params: ['ldsc_ld'],              ancestry_params: ['ldsc_ld'],  conda_ok: false ],
        lava: [ scope: 'local',  required_params: ['lava_ref', 'lava_loci'], ancestry_params: ['lava_ref'], conda_ok: false ],
        sumher: [ scope: 'global', required_params: ['sumher_tagfile'],     ancestry_params: ['sumher_tagfile'], conda_ok: true ],
        // Popcorn estimates the CROSS-ancestry pairs (the others take same-ancestry pairs)
        popcorn: [ scope: 'global', required_params: ['popcorn_cfile'],     ancestry_params: [],                 conda_ok: false, cross_ancestry: true ],
    ]
}

// GWAS ancestry of a trait (samplesheet column 'ancestry'; default EUR)
def traitAncestry(Map meta) {
    return (meta.ancestry ?: 'EUR').toString()
}

// Parameter holding a per-ancestry reference: ldsc_ld (EUR), ldsc_ld_eas (EAS), ...
def ancestryParam(String name, String ancestry) {
    return ancestry == 'EUR' ? name : "${name}_${ancestry.toLowerCase()}".toString()
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
    def before_text = gencorrLogo(monochrome_logs)
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
    // A GWAS Catalog accession in 'sumstats' becomes meta.gwas_catalog (the
    // file is downloaded by GWASCATALOG_FETCH) and sumstats [].
    //
    def samples = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, sumstats -> resolveSumstats(cleanMeta(meta), sumstats as String) }
    validateSamples(samples)

    //
    // Multivariate GWAMA: each group of --gwama_groups becomes a new trait, made
    // from its traits' harmonized files and analysed like a samplesheet trait.
    // [ [ group_meta, [ member metas ] ], ... ]
    //
    def sample_meta = samples.collectEntries { row -> [ (row[0].id): row[0] ] }
    def gwama_list = params.gwama_groups
        ? makeGwamaGroups(readGwamaGroups(params.gwama_groups as String), sample_meta)
        : []
    validateGwamaReferences(params, workflow.profile, gwama_list.collect { g -> traitAncestry(g[0]) }.unique().sort())

    //
    // Trait pairs, decided once here so every method analyses the same pairs.
    //
    def ids   = samples.collect { row -> row[0].id } + gwama_list.collect { g -> g[0].id }
    def pair_list = makePairs(ids, params.pair_mode, params.targets,
                              params.pairs ? readPairsFile(params.pairs) : null)
    def n_gwama = gwama_list ? " + ${gwama_list.size()} combined (GWAMA)" : ''
    log.info "gencorr: ${samples.size()} traits${n_gwama}, ${pair_list.size()} pairs (pair_mode = ${params.pair_mode}), methods: ${methods.join(', ')}"

    def meta_by_id = sample_meta + gwama_list.collectEntries { g -> [ (g[0].id): g[0] ] }

    //
    // Ancestry: LDSC, LAVA and SumHer analyse same-ancestry pairs with that
    // ancestry's references and report cross-ancestry pairs with status
    // cross_ancestry; Popcorn estimates the cross-ancestry pairs.
    //
    def same_ancestry = pair_list.findAll { p -> traitAncestry(meta_by_id[p[0]]) == traitAncestry(meta_by_id[p[1]]) }
    def n_cross       = pair_list.size() - same_ancestry.size()
    if (n_cross) {
        def by_popcorn = 'popcorn' in methods ? ' They are estimated by Popcorn.' : ' Add --methods popcorn (with --popcorn_cfile) to estimate them across ancestries.'
        log.warn "${n_cross} of ${pair_list.size()} pairs have traits of different ancestry; LDSC, LAVA and SumHer assume one population, so they report these pairs with status 'cross_ancestry'.${by_popcorn}"
    }
    validateReferences(params, methods, same_ancestry.collect { p -> traitAncestry(meta_by_id[p[0]]) }.unique().sort())
    if ('popcorn' in methods) {
        def cross = pair_list.findAll { p -> !(p in same_ancestry) }
        validatePopcornAncestries(params.popcorn_cfile_ancestries as String,
            cross.collect { p -> [ traitAncestry(meta_by_id[p[0]]), traitAncestry(meta_by_id[p[1]]) ] })
    }

    emit:
    samplesheet = channel.fromList(samples)                                        // [ meta, sumstats ]
    pairs       = channel.fromList(pair_list.collect { p -> [ meta_by_id[p[0]], meta_by_id[p[1]] ] }) // [ meta1, meta2 ]
    gwama       = channel.fromList(gwama_list)                                        // [ group_meta, [ member metas ] ]
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
// Start-up banner, printed before the parameter summary and the --help text.
// Plain ASCII so it renders in any terminal and in cluster log files.
//
def gencorrLogo(monochrome_logs = true) {
    def c   = logColours(monochrome_logs) as Map
    def art = [
        '   __ _  ___ _ __   ___ ___  _ __ _ __',
        '  / _` |/ _ \\ \'_ \\ / __/ _ \\| \'__| \'__|',
        ' | (_| |  __/ | | | (_| (_) | |  | |',
        '  \\__, |\\___|_| |_|\\___\\___/|_|  |_|',
    ]
    def lines = art.collect { line -> "${c.green}${line}${c.reset}".toString() }
    lines << "${c.green}  |___/${c.reset}   ${c.dim}global + local genetic correlation${c.reset}".toString()
    lines << "  ${c.bold}${workflow.manifest.name} ${workflow.manifest.version}${c.reset}".toString()
    return '\n' + lines.join('\n') + "\n${c.dim}${'-' * 60}${c.reset}\n"
}

//
// Drop null / empty-list / empty-string meta values produced by empty samplesheet cells
//
def cleanMeta(Map meta) {
    return meta.findAll { _k, v -> v != null && v != [] && v != '' }
}

//
// samplesheet 'sumstats': a GWAS Catalog accession (GCST + digits) or a file
// path / URL. Returns [ meta, file ] or [ meta + [gwas_catalog: accession], [] ].
//
def resolveSumstats(Map meta, String sumstats) {
    def value = sumstats.trim()
    if (value ==~ /GCST\d+/) {
        return [ meta + [ gwas_catalog: value ], [] ]
    }
    def f = file(value)
    if (!f.exists()) {
        error("Samplesheet row '${meta.id}': sumstats '${value}' is neither an existing file nor a GWAS Catalog accession (e.g. GCST90000016).")
    }
    return [ meta, f ]
}

//
// Samplesheet consistency checks that the JSON schema cannot express
//
def validateSamples(List samples) {
    samples.each { meta, _sumstats ->
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
//   target : each target trait against every other trait ("one vs group")
//   pairs  : exactly the pairs listed in a CSV (columns trait1,trait2)
// Pairs are unordered: [a, b] and [b, a] are the same pair; the first seen wins.
// (Written without switch/for loops: Nextflow >= 26 parses strict syntax.)
//
def makePairs(List ids, String mode, String targets, List listed) {
    def candidates = []
    if (mode == 'all') {
        candidates = (0..<ids.size()).collectMany { i ->
            ((i + 1)..<ids.size()).collect { j -> [ ids[i], ids[j] ] }
        }
    }
    else if (mode == 'target') {
        if (!targets) {
            error("--pair_mode target requires --targets (comma-separated trait IDs).")
        }
        def tlist   = targets.tokenize(',').collect { t -> t.trim() }
        def unknown = tlist.findAll { t -> !(t in ids) }
        if (unknown) {
            error("--targets: unknown trait ID(s) ${unknown}. Known: ${ids}")
        }
        candidates = tlist.collectMany { t -> ids.findAll { o -> o != t }.collect { o -> [ t, o ] } }
    }
    else if (mode == 'pairs') {
        if (listed == null) {
            error("--pair_mode pairs requires --pairs <csv with trait1,trait2>.")
        }
        listed.each { p ->
            def unknown = p.findAll { x -> !(x in ids) }
            if (unknown) {
                error("--pairs: unknown trait ID(s) ${unknown} in pair ${p}. Known: ${ids}")
            }
        }
        candidates = listed.collect { p -> [ p[0], p[1] ] }
    }
    else {
        error("--pair_mode must be one of all, target, pairs (got '${mode}').")
    }

    def self_pairs = candidates.findAll { p -> p[0] == p[1] }
    if (self_pairs) {
        error("Pair '${self_pairs[0][0]}' vs '${self_pairs[0][1]}': a trait cannot be paired with itself.")
    }
    // drop reversed duplicates, keep first occurrence and its orientation
    def out = candidates.unique(false) { p -> [ p[0], p[1] ].sort().join('\u0000') }
    if (!out) {
        error("No trait pairs to analyse with pair_mode '${mode}'.")
    }
    return out
}

def readPairsFile(String path) {
    def lines  = file(path, checkIfExists: true).readLines().findAll { line -> line.trim() }
    def header = lines[0].split(/[,\t]/).collect { x -> x.trim() }
    if (header.take(2) != ['trait1', 'trait2']) {
        error("--pairs: header must start with trait1,trait2 (got ${header})")
    }
    return lines.drop(1).collect { line -> line.split(/[,\t]/).collect { x -> x.trim() }.take(2) }
}

//
// --gwama_groups: CSV (or TSV) with header group,trait and an optional sign
// column. Returns the rows as [ group, trait, sign ] (sign '' when empty).
//
def readGwamaGroups(String path) {
    def lines  = file(path, checkIfExists: true).readLines().findAll { line -> line.trim() }
    def header = lines[0].split(/[,\t]/).collect { x -> x.trim() }
    if (header.take(2) != ['group', 'trait'] || (header.size() > 2 && header[2] != 'sign')) {
        error("--gwama_groups: header must be group,trait or group,trait,sign (got ${header})")
    }
    return lines.drop(1).collect { line ->
        def cells = line.split(/[,\t]/, -1).collect { x -> x.trim() }
        [ cells[0], cells.size() > 1 ? cells[1] : '', cells.size() > 2 ? cells[2] : '' ]
    }
}

//
// GWAMA groups from the rows of --gwama_groups ([ group, trait, sign ]) and the
// samplesheet metas by ID. Returns [ [ group_meta, [ member metas ] ], ... ] in
// file order; group_meta = [ id, ancestry, gwama_traits, gwama_signs ].
// A group needs >= 2 different samplesheet traits of one ancestry; its name
// follows the samplesheet ID rules and must not be a trait ID.
//
def makeGwamaGroups(List rows, Map meta_by_id) {
    def by_group = [:]
    rows.each { row ->
        def group = row[0]
        def trait = row[1]
        def sign  = row[2] ?: '1'
        if (!(group ==~ /^[A-Za-z](_?[A-Za-z0-9.-])*$/)) {
            error("--gwama_groups: group name '${group}' must start with a letter and contain only letters, digits, '.', '-' and single '_' (not at the end).")
        }
        if (meta_by_id.containsKey(group)) {
            error("--gwama_groups: group name '${group}' is also a trait ID in the samplesheet; the combined trait needs its own name.")
        }
        if (!meta_by_id.containsKey(trait)) {
            error("--gwama_groups: group '${group}': unknown trait '${trait}'. Known: ${meta_by_id.keySet() as List}")
        }
        if (!(sign in ['1', '+1', '-1'])) {
            error("--gwama_groups: group '${group}', trait '${trait}': sign must be 1 or -1 (got '${sign}').")
        }
        if (!by_group.containsKey(group)) {
            by_group[group] = []
        }
        if (by_group[group].any { m -> m[0] == trait }) {
            error("--gwama_groups: group '${group}' lists trait '${trait}' twice.")
        }
        by_group[group] << [ trait, sign == '-1' ? -1 : 1 ]
    }
    if (!by_group) {
        error("--gwama_groups: the file lists no groups.")
    }
    return by_group.collect { group, members ->
        if (members.size() < 2) {
            error("--gwama_groups: group '${group}' has ${members.size()} trait; N-GWAMA combines at least 2.")
        }
        def metas      = members.collect { m -> meta_by_id[m[0]] }
        def ancestries = metas.collect { m -> traitAncestry(m) }.unique()
        if (ancestries.size() > 1) {
            error("--gwama_groups: group '${group}' mixes ancestries ${ancestries}; N-GWAMA assumes one population (LDSC intercepts, LD).")
        }
        def gmeta = [ id: group, ancestry: ancestries[0],
                      gwama_traits: members.collect { m -> m[0] }, gwama_signs: members.collect { m -> m[1] } ]
        [ gmeta, metas ]
    }
}

//
// GWAMA takes h2 and intercepts from LDSC: every group's ancestry needs LD
// scores, and LDSC has no conda package.
//
def validateGwamaReferences(Map p, String profile, List ancestries) {
    if (!ancestries) {
        return
    }
    ancestries.each { a ->
        def name = ancestryParam('ldsc_ld', a)
        if (!p[name]) {
            error("--gwama_groups: combining ${a} traits needs LDSC LD scores (--${name}) for their h2 and intercepts.")
        }
    }
    if ((profile ?: '').tokenize(',').intersect(['conda', 'mamba'])) {
        error("--gwama_groups runs LDSC, which has no conda package. Use -profile docker, singularity or apptainer.")
    }
}

//
// Cross-parameter checks. Returns the list of methods to run.
//
def validateGencorrParams(Map p, String profile) {
    def registry = methodRegistry()
    def methods  = (p.methods ?: '').tokenize(',').collect { m -> m.trim() }.unique()
    def unknown  = methods.findAll { m -> !registry.containsKey(m) }
    if (!methods) {
        error("--methods is empty. Choose from: ${registry.keySet().join(', ')}")
    }
    if (unknown) {
        error("--methods: unknown method(s) ${unknown}. Choose from: ${registry.keySet().join(', ')}")
    }

    methods.each { m ->
        // per-ancestry references are checked later, once the samplesheet is read
        registry[m].required_params.findAll { rp -> !(rp in registry[m].ancestry_params) }.each { rp ->
            if (!p[rp]) {
                error("--methods ${m} requires --${rp}.")
            }
        }
        def conda = (profile ?: '').tokenize(',').intersect(['conda', 'mamba'])
        if (conda && !registry[m].conda_ok) {
            error("Method '${m}' has no conda package (it is installed from GitHub in its container). Use -profile docker, singularity or apptainer.")
        }
    }
    if ('lava' in methods && p.lava_overlap == 'ldsc' && !('ldsc' in methods)) {
        error("--lava_overlap ldsc estimates sample overlap from LDSC intercepts, so it needs --methods to include ldsc. Alternatively use --lava_overlap none or a matrix file.")
    }
    if (!(p.lava_snp_set in ['pairwise', 'joint'])) {
        error("--lava_snp_set must be 'pairwise' or 'joint' (got '${p.lava_snp_set}').")
    }
    if ((p.ldsc_ld || p.ldsc_ld_eas) && !p.ldsc_merge_alleles && 'ldsc' in methods) {
        log.warn "No --ldsc_merge_alleles given. LDSC fails a pair when an rsID has different alleles in the two traits; the usual choice is w_hm3.snplist."
    }
    return methods
}

//
// Per-ancestry references: every ancestry with at least one same-ancestry pair
// needs its reference for every selected method (e.g. EAS + ldsc -> --ldsc_ld_eas).
//
def validateReferences(Map p, List methods, List ancestries) {
    def registry = methodRegistry()
    methods.each { m ->
        registry[m].ancestry_params.each { rp ->
            ancestries.each { a ->
                def name = ancestryParam(rp, a)
                if (!p[name]) {
                    error("--methods ${m}: traits with ancestry ${a} are paired with each other, which needs --${name}.")
                }
            }
        }
    }
}

//
// Popcorn: --popcorn_cfile_ancestries names the score file's populations 1 and 2.
// It must be two different ancestries, and every cross-ancestry pair must be
// one trait of each (METHOD_POPCORN orders the pair by it; a wrong value would
// give Popcorn the traits in the wrong populations without any error).
// pair_ancestries: [ [ancestry1, ancestry2], ... ] of the cross-ancestry pairs.
// Returns the two ancestries.
//
def validatePopcornAncestries(String cfile_ancestries, List pair_ancestries) {
    def pops = (cfile_ancestries ?: '').tokenize(',').collect { a -> a.trim() }
    if (pops.size() != 2 || pops[0] == pops[1]) {
        error("--popcorn_cfile_ancestries must name two different ancestries, in the score file's order, e.g. 'EUR,EAS' (got '${cfile_ancestries}').")
    }
    def other = pair_ancestries.find { pair -> (pair as Set) != (pops as Set) }
    if (other) {
        error("--methods popcorn: a ${other[0]} x ${other[1]} pair cannot be estimated with a score file for ${pops[0]} and ${pops[1]} (--popcorn_cfile_ancestries).")
    }
    return pops
}
