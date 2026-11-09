#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Regenerates the small LDSC LD-score reference committed under tests/data/ldsc/.
#
# Provenance
#   Genotypes: LAVA vignette reference (1000 Genomes EUR subset, GRCh37),
#              josefin-werme/LAVA @ e729a245f7b6923967a96804fbf5246eadf2d6c6
#              vignettes/data/g1000_test.{bed,bim,fam}
#   Software : plink2 (bioconda), ldsc 1.0.1 (bioconda, python 2.7)
#
# Only the derived LD scores are committed. The genotype and summary statistic
# files are fetched at test time straight from the pinned LAVA commit
# (see conf/test.config), so nothing from the LAVA repository is redistributed.
#
# These LD scores cover ~350k SNPs in 42 small regions. They are only for
# exercising the pipeline. Never use them for real analyses: use the
# published LD scores for your ancestry (e.g. eur_w_ld_chr, UKBB/Pan-UKB).
#
# Usage: bash tests/data/make_test_data.sh <ldsc_env_bin> <plink2_bin>
# -----------------------------------------------------------------------------
set -euo pipefail

LDSC_BIN=${1:?path to directory containing ldsc.py}
PLINK2=${2:?path to plink2 binary}
COMMIT=e729a245f7b6923967a96804fbf5246eadf2d6c6
BASE=https://raw.githubusercontent.com/josefin-werme/LAVA/${COMMIT}/vignettes/data
HERE=$(cd "$(dirname "$0")" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

for ext in bed bim fam; do
    curl -sSfL -o "$TMP/g1000_test.$ext" "$BASE/g1000_test.$ext"
done

mkdir -p "$HERE/ldsc"
for chr in $(cut -f1 "$TMP/g1000_test.bim" | sort -un); do
    "$PLINK2" --bfile "$TMP/g1000_test" --chr "$chr" --make-bed \
        --out "$TMP/chr$chr" --silent
    # 1 Mb window (the test bim has no cM positions). --yes-really is needed because
    # each chromosome here is a few small regions, which ldsc mistakes for a
    # whole-chromosome window; it is cheap at this size.
    "$LDSC_BIN/python" "$LDSC_BIN/ldsc.py" --bfile "$TMP/chr$chr" --l2 \
        --ld-wind-kb 1000 --yes-really --out "$HERE/ldsc/$chr" > /dev/null
    rm -f "$HERE/ldsc/$chr.log"
done
# The test data covers 14 chromosomes, but --ref-ld-chr expects all of 1-22,
# so merge everything into one LD-score set used with --ref-ld (no -chr).
# This also exercises the pipeline's ldsc_ld_split_by_chr = false code path.
{
    zcat "$HERE"/ldsc/[0-9]*.l2.ldscore.gz | head -1
    for f in $(ls "$HERE"/ldsc/[0-9]*.l2.ldscore.gz | sort -V); do zcat "$f" | tail -n +2; done
} | gzip -n > "$HERE/ldsc/all.l2.ldscore.gz"
for m in M M_5_50; do
    cat "$HERE"/ldsc/[0-9]*.l2.$m | awk '{s+=$1} END{print s}' > "$HERE/ldsc/all.l2.$m"
done
rm -f "$HERE"/ldsc/[0-9]*.l2.*
echo "LD scores written to $HERE/ldsc (prefix: all)"
