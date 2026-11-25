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
        // (score files per ancestry pair are checked in popcornScoreSources)
        popcorn: [ scope: 'global', required_params: [],                    ancestry_params: [],                 conda_ok: false, cross_ancestry: true ],
    ]
}

// GWAS ancestry of a trait (samplesheet column 'ancestry'; default EUR)
def traitAncestry(Map meta) {
    return (meta.ancestry ?: 'EUR').toString()
}

// Ancestries with their own references (one parameter set each, see ancestryParam).
// The samplesheet also accepts MULTI: a multi-ancestry meta-analysis has no
// matching LD reference, so it is harmonized (QC report) but never analysed.
def supportedAncestries() {
    return ['EUR', 'EAS', 'AFR', 'SAS', 'AMR']
}

// Pairs with a MULTI trait get status 'multi_ancestry' from every method
def isMultiAncestryPair(Map m1, Map m2) {
    return traitAncestry(m1) == 'MULTI' || traitAncestry(m2) == 'MULTI'
}

// Parameter holding a per-ancestry reference: ldsc_ld (EUR), ldsc_ld_eas (EAS), ...
def ancestryParam(String name, String ancestry) {
    return ancestry == 'EUR' ? name : "${name}_${ancestry.toLowerCase()}".toString()
}

// Per-ancestry setting that falls back to the EUR one when not set
// (e.g. ldsc_ld_prefix_afr -> ldsc_ld_prefix)
def ancestrySetting(Map p, String name, String ancestry) {
    def value = p[ancestryParam(name, ancestry)]
    return value != null ? value : p[name]
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
    // Trait pairs, decided once here so every method analyses the same pairs.
    //
    def ids   = samples.collect { row -> row[0].id }
    def pair_list = makePairs(ids, params.pair_mode, params.targets,
                              params.pairs ? readPairsFile(params.pairs) : null)
    log.info "gencorr: ${samples.size()} traits, ${pair_list.size()} pairs (pair_mode = ${params.pair_mode}), methods: ${methods.join(', ')}"

    def meta_by_id = samples.collectEntries { row -> [ (row[0].id): row[0] ] }

    //
    // Ancestry: LDSC, LAVA and SumHer analyse same-ancestry pairs with that
    // ancestry's references and report cross-ancestry pairs with status
    // cross_ancestry; Popcorn estimates the cross-ancestry pairs.
    //
    def multi_pairs   = pair_list.findAll { p -> isMultiAncestryPair(meta_by_id[p[0]], meta_by_id[p[1]]) }
    if (multi_pairs) {
        log.warn "${multi_pairs.size()} of ${pair_list.size()} pairs involve a trait with ancestry MULTI (multi-ancestry meta-analysis). No LD reference matches such a sample, so every method reports these pairs with status 'multi_ancestry'. Use the ancestry-specific GWAS instead."
    }
    def analysable    = pair_list.findAll { p -> !(p in multi_pairs) }
    def same_ancestry = analysable.findAll { p -> traitAncestry(meta_by_id[p[0]]) == traitAncestry(meta_by_id[p[1]]) }
    def n_cross       = analysable.size() - same_ancestry.size()
    if (same_ancestry.any { p -> traitAncestry(meta_by_id[p[0]]) == 'AMR' }) {
        log.warn "Some same-ancestry pairs are AMR (admixed). LD scores from an external reference underestimate h2 in admixed samples; prefer covariate-adjusted LD scores (Pan-UKB AMR) and treat AMR results as exploratory (docs/usage.md)."
    }
    if (n_cross) {
        def by_popcorn = 'popcorn' in methods ? ' They are estimated by Popcorn.' : ' Add --methods popcorn to estimate them across ancestries.'
        log.warn "${n_cross} of ${pair_list.size()} pairs have traits of different ancestry; LDSC, LAVA and SumHer assume one population, so they report these pairs with status 'cross_ancestry'.${by_popcorn}"
    }
    validateReferences(params, methods, same_ancestry.collect { p -> traitAncestry(meta_by_id[p[0]]) }.unique().sort())
    if ('popcorn' in methods) {
        def cross_ancestries = analysable.findAll { p -> !(p in same_ancestry) }
            .collect { p -> [ traitAncestry(meta_by_id[p[0]]), traitAncestry(meta_by_id[p[1]]) ] }
        if (cross_ancestries.any { pair -> 'AMR' in pair }) {
            log.warn "Popcorn assumes non-admixed samples, so cross-ancestry pairs with an AMR trait are reported with status 'admixed_ancestry'."
        }
        // stops the run when an ancestry pair has no score file and no references
        popcornScoreSources(params, params.popcorn_cfiles ? readPopcornCfiles(params.popcorn_cfiles as String) : [], cross_ancestries)
    }

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
    if (supportedAncestries().any { a -> p[ancestryParam('ldsc_ld', a)] } && !p.ldsc_merge_alleles && 'ldsc' in methods) {
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
                    def hint = (name == 'sumher_tagfile_amr')
                        ? ' LDAK publishes no AMR tagging file: compute one with ldak --calc-tagging from an AMR reference, or leave sumher out of --methods.'
                        : ''
                    error("--methods ${m}: traits with ancestry ${a} are paired with each other, which needs --${name}.${hint}")
                }
            }
        }
    }
}

//
// Popcorn: --popcorn_cfile_ancestries names the populations 1 and 2 of
// --popcorn_cfile, in the score file's order (METHOD_POPCORN orders each pair by
// it; a wrong value would give Popcorn the traits in the wrong populations
// without any error). Returns the two ancestries.
//
def validatePopcornAncestries(String cfile_ancestries) {
    def pops = (cfile_ancestries ?: '').tokenize(',').collect { a -> a.trim() }
    if (pops.size() != 2 || pops[0] == pops[1]) {
        error("--popcorn_cfile_ancestries must name two different ancestries, in the score file's order, e.g. 'EUR,EAS' (got '${cfile_ancestries}').")
    }
    return pops
}

//
// --popcorn_cfiles: CSV with header ancestry_1,ancestry_2,cfile (one pre-computed
// score file per ancestry pair, populations in the file's order).
// Returns [ [ancestry_1, ancestry_2, cfile], ... ].
//
def readPopcornCfiles(String path) {
    def lines  = file(path, checkIfExists: true).readLines().findAll { line -> line.trim() }
    def header = lines ? lines[0].split(',').collect { x -> x.trim() } : []
    if (header != ['ancestry_1', 'ancestry_2', 'cfile']) {
        error("--popcorn_cfiles: header must be ancestry_1,ancestry_2,cfile (got ${header})")
    }
    return lines.drop(1).collect { line ->
        def row = line.split(',').collect { x -> x.trim() }
        if (row.size() != 3 || !(row[0] in supportedAncestries()) || !(row[1] in supportedAncestries()) || row[0] == row[1]) {
            error("--popcorn_cfiles: '${line}' must be two different ancestries (${supportedAncestries().join(', ')}) and a file.")
        }
        row
    }
}

//
// Popcorn score file for each ancestry pair of the cross-ancestry trait pairs,
// first match wins:
//   1. a row of --popcorn_cfiles (cfile_rows) for the two ancestries
//   2. --popcorn_cfile when --popcorn_cfile_ancestries names the two ancestries
//   3. computed by POPCORN_COMPUTE from --popcorn_ref[_<x>] of both ancestries
//      (population 1 = the ancestry listed first in supportedAncestries())
// Pairs with an AMR trait are skipped: Popcorn assumes non-admixed samples, so
// they get status admixed_ancestry. A pair without any source stops the run.
// pair_ancestries: [ [a, b], ... ]. Returns one entry per distinct ancestry pair:
// [ ancestries: [pop1, pop2], cfile: <path> or null (= compute) ].
//
def popcornScoreSources(Map p, List cfile_rows, List pair_ancestries) {
    def single = p.popcorn_cfile ? validatePopcornAncestries(p.popcorn_cfile_ancestries as String) : null
    def order  = supportedAncestries()
    def wanted = pair_ancestries
        .findAll { pair -> !('AMR' in pair) }
        .collect { pair -> pair.collect { a -> a as String }.sort { a, b -> order.indexOf(a) <=> order.indexOf(b) } }
        .unique()
    return wanted.collect { pair ->
        def row = cfile_rows.find { r -> ([r[0], r[1]] as Set) == (pair as Set) }
        if (row) {
            return [ ancestries: [ row[0], row[1] ], cfile: row[2] ]
        }
        if (single && (single as Set) == (pair as Set)) {
            return [ ancestries: single, cfile: p.popcorn_cfile ]
        }
        def missing = pair.collect { a -> ancestryParam('popcorn_ref', a) }.findAll { name -> !p[name] }
        if (missing) {
            error("--methods popcorn: ${pair[0]} x ${pair[1]} pairs need a score file: a row in --popcorn_cfiles, --popcorn_cfile with --popcorn_cfile_ancestries ${pair.join(',')}, or PLINK references to compute one (${missing.collect { name -> '--' + name }.join(', ')}).")
        }
        return [ ancestries: pair, cfile: null ]
    }
}
