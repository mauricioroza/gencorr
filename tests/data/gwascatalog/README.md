# Synthetic GWAS Catalog mirror

A folder laid out like the GWAS Catalog FTP
(`<range>/<accession>/harmonised/<accession>.h.tsv.gz` and `-meta.yaml`), used
by the `GWASCATALOG_FETCH` module test with `--base tests/data/gwascatalog`.
GCST90000999 here is **not** the real study of that accession: three made-up
variants and a made-up sample size. Regenerate with the recipe in
`tests/data/make_test_data.sh`.
