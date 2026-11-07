# Adding a method

Every method plugs into the same two contracts, so the rest of the pipeline
(harmonization, pairing, merging, comparison tool) works for it unchanged:

| contract | defined in | what it fixes |
|---|---|---|
| harmonized input | `bin/gencorr_common.py` → `HARMONIZED_COLUMNS` | SNP, CHR, BP, A1 (effect), A2, Z, BETA, SE, P, N, EAF, INFO |
| standard output | `bin/gencorr_common.py` → `RESULT_COLUMNS` | one row per (pair, locus); method-specific values go in the JSON `extra` column |

A method subworkflow therefore always looks like this:

```nextflow
workflow METHOD_X {
    take:
    ch_harmonized   // [ meta, harmonized.tsv.gz ]   every trait
    ch_pairs        // [ meta1, meta2 ]              the pairs to analyse
    ...             // the method's own reference data
    main:
    EXPORT_SUMSTATS ( ch_harmonized, 'x' )           // 1. format for the tool
    X_RUN ( ... )                                    // 2. run the tool
    X_PARSE ( ... )                                  // 3. -> standard results table
    emit:
    results = X_PARSE.out.results                    // path: standard table
}
```

## Worked example: SumHer (LDAK) global rg

BigaGWAS already runs SumHer, so this is the natural next one. Five steps:

### 1. Input format: one entry in `bin/export_sumstats.py`

SumHer's `--summary` file wants `Predictor A1 A2 n Z`:

```python
FORMATS["sumher"] = {
    "columns": [("Predictor", "SNP"), ("A1", "A1"), ("A2", "A2"), ("n", "N"), ("Z", "Z")],
    "optional": [],
    "min_coverage": 1.0,
    "sep": "\t",
    "compress": False,
}
```

### 2. Modules: `modules/local/sumher/rg/main.nf` (+ `environment.yml`)

Copy `modules/local/ldsc/rg/main.nf` as a template. Keep:
`tag`, `label`, `conda`, `container` (pin an exact version; LDAK 6 is on
bioconda as `ldak6`), the `versions` topic line, `ext.args`, and a `stub:`
block. The command is the one from BigaGWAS'
`static/shell/specific_sumher.sh`:

```bash
ldak --sum-cors ${meta.id} --summary ${s1} --summary2 ${s2} \
     --tagfile ${tagfile} --allow-ambiguous YES --check-sums NO ${args}
```

### 3. Parser: `bin/parse_sumher.py`

Read the `.cors` output and write `RESULT_COLUMNS` with
`gencorr_common.validate_results()` + `write_tsv()`. Copy the structure of
`bin/parse_ldsc.py`: take every number from the tool's own output, put
anything that has no standard column into `extra`, and give every pair a
`status`. Add unit tests with a real output file as fixture in
`tests/python/fixtures/sumher/`.

### 4. Subworkflow: `subworkflows/local/method_sumher/main.nf`

Same `take:`/`emit:` as `METHOD_LDSC`. Reuse the pair-batching code from
`METHOD_LDSC` if the tool accepts several traits per call.

### 5. Register it (two places)

`subworkflows/local/utils_nfcore_gencorr_pipeline/main.nf`:

```groovy
def methodRegistry() {
    return [
        ldsc:   [ scope: 'global', required_params: ['ldsc_ld'],               conda_ok: true  ],
        lava:   [ scope: 'local',  required_params: ['lava_ref', 'lava_loci'], conda_ok: false ],
        sumher: [ scope: 'global', required_params: ['sumher_tagfile'],        conda_ok: true  ],
    ]
}
```

`workflows/gencorr.nf`:

```groovy
if ('sumher' in methods) {
    METHOD_SUMHER ( ch_harmonized, ch_pairs, channel.value(file(params.sumher_tagfile)) )
    ch_results = ch_results.mix( METHOD_SUMHER.out.results )
}
```

Then add `sumher_*` params to `nextflow.config`, describe them in
`nextflow_schema.json` (`nf-core pipelines schema build`), and add publishing
rules to `conf/modules.config`.

## Checklist before merging a new method

- [ ] Module(s) with pinned container + conda env, `versions` topic, `stub:`
- [ ] Parser writes `RESULT_COLUMNS` only (`validate_results` passes)
- [ ] Module test next to each new module (`modules/local/<name>/tests/main.nf.test`,
      real-data test + `-stub` test, tags `modules`, `modules_local`, `<name>`);
      see `modules/local/ldsc/rg/tests/` for the `setup` pattern
- [ ] Unit tests for the parser on real tool output (`tests/python/`)
- [ ] **Oracle test**: the tool run by hand on test data gives the same numbers
      as the pipeline (see `tests/python/test_ldsc_oracle.py` for the pattern)
- [ ] `conf/test.config` exercises it; `nf-test test tests/default.nf.test`
      passes and snapshots are updated (`--update-snapshot`) with a reviewed diff
- [ ] Docs: parameter table in `usage.md`, columns/`extra` keys in `output.md`,
      citation in `CITATIONS.md`

## Candidate methods

| method | scope | tool | note |
|---|---|---|---|
| SumHer | global | LDAK | in BigaGWAS; tagging files per ancestry |
| Popcorn | global, cross-ancestry | Popcorn | in BigaGWAS |
| HDL | global | R package | lower variance than LDSC |
| HDL-L | local | R package | random-effects local rg, a useful cross-check of LAVA |
| SUPERGNOVA, rho-HESS | local | | alternatives with different models |
