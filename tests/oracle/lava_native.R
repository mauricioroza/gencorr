#!/usr/bin/env Rscript
# =============================================================================
# Oracle for tests: runs LAVA the way its own vignette does, with NO gencorr
# code involved: process.input() -> process.locus() -> run.univ.bivar().
#
#   lava_native.R <info> <ref_prefix> <loci> <out.tsv> <seed> <univ_thresh> <pheno1> <pheno2>
#
# The only addition is set.seed() before each locus, using the same seed
# formula as bin/run_lava.R, so Monte Carlo CIs/p-values are comparable
# draw-for-draw. If gencorr reproduces this table exactly, the runner adds
# nothing to (and removes nothing from) LAVA's own computation.
# =============================================================================
suppressPackageStartupMessages(library(LAVA))
a <- commandArgs(trailingOnly = TRUE)
info <- a[1]; ref <- a[2]; loci_file <- a[3]; out <- a[4]
seed <- as.numeric(a[5]); thresh <- as.numeric(a[6]); phenos <- a[7:8]

seed_for <- function(base, ...) {   # identical to bin/run_lava.R
    h <- 0
    for (code in utf8ToInt(paste(..., sep = "|"))) h <- (h * 31 + code) %% 2147483647
    as.integer((base + h) %% 2147483647)
}

input <- process.input(input.info.file = info, sample.overlap.file = NULL,
                       ref.prefix = ref, phenos = phenos)
loci <- read.loci(loci_file)
rows <- list()
for (i in seq_len(nrow(loci))) {
    locus <- process.locus(loci[i, ], input)
    if (is.null(locus)) next
    set.seed(seed_for(seed, loci$LOC[i], phenos[1], phenos[2]))
    res <- run.univ.bivar(locus, univ.thresh = thresh)
    u <- res$univ
    row <- data.frame(locus = as.character(loci$LOC[i]),
                      h2_1 = u$h2.obs[u$phen == phenos[1]], h2_2 = u$h2.obs[u$phen == phenos[2]],
                      univ_p_1 = u$p[u$phen == phenos[1]], univ_p_2 = u$p[u$phen == phenos[2]],
                      rg = NA, rg_ci_low = NA, rg_ci_high = NA, p = NA)
    if (!is.null(res$bivar)) {
        b <- res$bivar
        row$rg <- b$rho; row$rg_ci_low <- b$rho.lower; row$rg_ci_high <- b$rho.upper; row$p <- b$p
    }
    rows[[length(rows) + 1]] <- row
}
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
