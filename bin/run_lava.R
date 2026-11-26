#!/usr/bin/env Rscript
# =============================================================================
# run_lava.R - run LAVA univariate + bivariate local genetic correlation for a
# set of loci and a set of trait pairs, writing the gencorr standard results.
#
# Usage (all arguments are --key value):
#   run_lava.R --info input.info.txt --ref <reference prefix> --loci loci.txt
#              --pairs pairs.tsv --out <prefix>
#              [--overlap sample_overlap.txt] [--univ-thresh 0.05]
#              [--adap-thresh 1e-4,1e-6] [--param-lim 1.25]
#              [--p-values TRUE] [--cis TRUE] [--seed 20201231] [--snp-set all]
#              [--mhc-region 6:25000000-34000000] [--exclude-mhc FALSE]
#
# Inputs
#   --info     LAVA input info file: phenotype, cases, controls, prevalence, filename
#   --ref      reference prefix (PLINK .bed/.bim/.fam, or LAVA LD-format prefix)
#   --loci     LAVA locus file (LOC CHR START STOP); any subset/chunk of loci
#   --pairs    TSV with header trait1<TAB>trait2; the pairs to analyse
#   --overlap  sample-overlap (sampling correlation) matrix, or omit / "none"
#
# Behaviour (kept identical to LAVA's own run.univ.bivar(), see note 1):
#   * run.univ() on every locus; the bivariate test for a pair is only run when
#     BOTH traits have univariate p < --univ-thresh (LAVA's recommended filter).
#   * Pairs that are not tested still get a row, with status explaining why, so
#     a missing estimate is never silent.
#   * Every row records extra.in_mhc: the locus overlaps --mhc-region (extended
#     MHC, in the coordinates of the locus file; GRCh37 for LAVA's own file).
#     Long-range LD there inflates local estimates. With --exclude-mhc TRUE those
#     loci are not analysed (status excluded_mhc).
#
# Reproducibility (note 2):
#   LAVA's CIs and bivariate p-values are Monte Carlo (matrixsampling::rwishart).
#   Before each bivariate call the RNG is seeded from (--seed, locus ID, trait1,
#   trait2). The result for a (locus, pair) therefore does not depend on how
#   loci are chunked across jobs, on job order, or on which other pairs are run.
#
# Notes
#   1. run.bivar() is called with exactly two phenotypes. Called with more, LAVA
#      draws the CIs for all pairs from one joint simulation, which would tie a
#      pair's CI to the other traits present.
#   2. Within one task, every trait listed in --info is loaded with
#      process.input(), and LAVA keeps only SNPs present in ALL loaded traits.
#      The Nextflow subworkflow controls this: lava_snp_set = 'pairwise' loads
#      2 traits per task, 'joint' loads every trait.
#   3. Loci without any analysis SNP are skipped BEFORE process.locus() and get
#      status "no_snps_in_locus". process.locus() scans the whole reference SNP
#      table (twice) for every locus, ~1 s per locus with the full 1000G EUR
#      reference, only to return NULL for such loci; on sparse or regional input
#      most of the 2,495 loci are empty. The count uses the same rule as LAVA
#      (reference CHR/POS within START..STOP, intersected with analysis SNPs),
#      so no locus that LAVA could analyse is skipped.
# =============================================================================

suppressPackageStartupMessages(library(LAVA))

# ---------------------------------------------------------------- arguments --
parse_args <- function(argv) {
    if (length(argv) %% 2 != 0) stop("arguments must be --key value pairs")
    keys <- sub("^--", "", argv[c(TRUE, FALSE)])
    vals <- argv[c(FALSE, TRUE)]
    as.list(setNames(vals, keys))
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(key, default = NULL, required = FALSE) {
    if (!is.null(opt[[key]])) return(opt[[key]])
    if (required) stop("missing required argument --", key)
    default
}
as_flag <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES")

info_file   <- get_opt("info", required = TRUE)
ref_prefix  <- get_opt("ref", required = TRUE)
loci_file   <- get_opt("loci", required = TRUE)
pairs_file  <- get_opt("pairs", required = TRUE)
out_prefix  <- get_opt("out", required = TRUE)
overlap     <- get_opt("overlap", "none")
univ_thresh <- as.numeric(get_opt("univ-thresh", "0.05"))
adap_thresh <- as.numeric(strsplit(get_opt("adap-thresh", "1e-4,1e-6"), ",")[[1]])
param_lim   <- as.numeric(get_opt("param-lim", "1.25"))
do_p        <- as_flag(get_opt("p-values", "TRUE"))
do_ci       <- as_flag(get_opt("cis", "TRUE"))
base_seed   <- as.numeric(get_opt("seed", "20201231"))
snp_set     <- get_opt("snp-set", "all")   # recorded only: "all" or the SNP-list file name
mhc_region  <- get_opt("mhc-region", "6:25000000-34000000")
exclude_mhc <- as_flag(get_opt("exclude-mhc", "FALSE"))
mhc <- regmatches(mhc_region, regexec("^(?:chr)?([0-9XYxy]+):([0-9]+)-([0-9]+)$", mhc_region, perl = TRUE))[[1]]
if (length(mhc) != 4) stop("--mhc-region must look like 6:25000000-34000000, got ", mhc_region)
mhc <- list(chr = toupper(mhc[2]), start = as.numeric(mhc[3]), stop = as.numeric(mhc[4]))
in_mhc <- function(loc) {
    !is.null(loc$CHR) && !is.null(loc$START) && !is.null(loc$STOP) &&
        toupper(sub("^chr", "", as.character(loc$CHR))) == mhc$chr &&
        loc$STOP >= mhc$start && loc$START <= mhc$stop
}
if (overlap %in% c("none", "NA", "")) overlap <- NULL

# ------------------------------------------------------------------ helpers --
# Deterministic, platform-independent seed from a base seed and strings.
# (Polynomial rolling hash mod 2^31-1 on UTF-8 code points.)
seed_for <- function(base, ...) {
    h <- 0
    for (code in utf8ToInt(paste(..., sep = "|"))) h <- (h * 31 + code) %% 2147483647
    as.integer((base + h) %% 2147483647)
}

json_value <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("null")
    if (is.numeric(x) && !is.finite(x)) return("null")   # Inf is not valid JSON
    if (is.numeric(x)) return(format(x, digits = 10, scientific = TRUE, trim = TRUE))
    if (is.logical(x)) return(tolower(as.character(x)))
    paste0("\"", gsub("[\"\\\\]", "", as.character(x)), "\"")
}
to_json <- function(lst) {
    paste0("{", paste0("\"", names(lst), "\":", vapply(lst, json_value, ""), collapse = ","), "}")
}

RESULT_COLUMNS <- c("method", "scope", "trait1", "trait2", "locus", "chr", "start",
                    "stop", "n_snps", "rg", "rg_se", "rg_ci_low", "rg_ci_high", "z",
                    "p", "h2_1", "h2_1_se", "h2_2", "h2_2_se", "status", "extra")

result_row <- function(loc, t1, t2, status, n_snps = NA, k = NA, u1 = NULL, u2 = NULL,
                       bv = NULL) {
    pick <- function(df, col) if (!is.null(df) && col %in% names(df)) df[[col]][1] else NA
    extra <- list(
        K = k,
        univ_p_1 = pick(u1, "p"), univ_p_2 = pick(u2, "p"),
        h2_latent_1 = pick(u1, "h2.latent"), h2_latent_2 = pick(u2, "h2.latent"),
        r2 = pick(bv, "r2"), r2_ci_low = pick(bv, "r2.lower"), r2_ci_high = pick(bv, "r2.upper"),
        univ_thresh = univ_thresh, seed = base_seed, snp_set = snp_set,
        in_mhc = in_mhc(loc)
    )
    data.frame(
        method = "lava", scope = "local", trait1 = t1, trait2 = t2,
        locus = as.character(loc$LOC), chr = as.character(loc$CHR),
        start = loc$START, stop = loc$STOP, n_snps = n_snps,
        rg = pick(bv, "rho"), rg_se = NA,
        rg_ci_low = pick(bv, "rho.lower"), rg_ci_high = pick(bv, "rho.upper"),
        z = NA, p = pick(bv, "p"),
        h2_1 = pick(u1, "h2.obs"), h2_1_se = NA, h2_2 = pick(u2, "h2.obs"), h2_2_se = NA,
        status = status, extra = to_json(extra), stringsAsFactors = FALSE
    )[, RESULT_COLUMNS]
}

# --------------------------------------------------------------------- main --
pairs <- read.delim(pairs_file, colClasses = "character")
stopifnot(all(c("trait1", "trait2") %in% names(pairs)), nrow(pairs) > 0)
phenos <- unique(c(pairs$trait1, pairs$trait2))

cat("LAVA version:", as.character(packageVersion("LAVA")), "\n")
cat("Traits loaded together (SNPs intersected across all of them):",
    paste(phenos, collapse = ", "), "\n")

# LAVA 0.1.5 uses the sample-overlap matrix without checking it. A block that is
# not positive definite (e.g. |r| >= 1 from LDSC intercepts) gives meaningless
# estimates, so stop with a clear message instead.
check_overlap <- function(path, phenos) {
    m <- as.matrix(read.table(path, header = TRUE, check.names = FALSE))
    use <- intersect(phenos, rownames(m))
    if (length(use) < 2) return(invisible(NULL))
    block <- m[use, use, drop = FALSE]
    min_eig <- min(eigen((block + t(block)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig <= 1e-8) {
        stop(sprintf(paste0(
            "sample-overlap matrix for traits %s is not positive definite (smallest eigenvalue %.3g). ",
            "LAVA would give meaningless estimates. Check the LDSC h2/gcov intercepts of these ",
            "traits (ldsc/gencorr.sample_overlap.txt), or run without --lava_overlap."),
            paste(use, collapse = ", "), min_eig), call. = FALSE)
    }
}
if (!is.null(overlap)) check_overlap(overlap, phenos)

input <- process.input(input.info.file = info_file, sample.overlap.file = overlap,
                       ref.prefix = ref_prefix, phenos = phenos)
loci <- read.loci(loci_file)
cat("Loci in this chunk:", nrow(loci), " Pairs:", nrow(pairs), "\n")

# ---- pre-count analysis SNPs per locus (note 3) ----------------------------- #
# Same definition as process.locus(): reference SNPs with CHR == locus CHR and
# START <= POS <= STOP, restricted to input$analysis.snps.
n_pre <- rep(NA_integer_, nrow(loci))
if (all(c("CHR", "START", "STOP") %in% names(loci)) && !("SNPS" %in% names(loci))) {
    ref_info <- input$reference$snp.info
    an <- ref_info[ref_info$SNP %in% input$analysis.snps, c("CHR", "POS")]
    pos_by_chr <- lapply(split(an$POS, as.character(an$CHR)), sort)
    n_pre <- vapply(seq_len(nrow(loci)), function(i) {
        pos <- pos_by_chr[[as.character(loci$CHR[i])]]
        if (is.null(pos)) return(0L)
        as.integer(findInterval(loci$STOP[i], pos) - findInterval(loci$START[i] - 1, pos))
    }, integer(1))
    cat("Loci with analysis SNPs:", sum(n_pre > 0), "of", nrow(loci),
        "(the others are skipped as no_snps_in_locus)\n")
}
t_start <- Sys.time()

results <- list()
univ_out <- list()
skip_mhc <- exclude_mhc & vapply(seq_len(nrow(loci)), function(i) in_mhc(loci[i, ]), logical(1))
if (any(skip_mhc)) cat("Loci in the MHC region", mhc_region, "excluded:", sum(skip_mhc), "
")
n_todo <- sum((is.na(n_pre) | n_pre > 0) & !skip_mhc)
n_done <- 0
for (i in seq_len(nrow(loci))) {
    loc <- loci[i, ]
    if (exclude_mhc && in_mhc(loc)) {
        for (j in seq_len(nrow(pairs)))
            results[[length(results) + 1]] <- result_row(loc, pairs$trait1[j], pairs$trait2[j],
                                                         "excluded_mhc")
        next
    }
    if (!is.na(n_pre[i]) && n_pre[i] == 0) {
        for (j in seq_len(nrow(pairs)))
            results[[length(results) + 1]] <- result_row(loc, pairs$trait1[j], pairs$trait2[j],
                                                         "no_snps_in_locus", n_snps = 0)
        next
    }
    n_done <- n_done + 1
    t_locus <- Sys.time()
    locus <- tryCatch(process.locus(loc, input), error = function(e) {
        message("locus ", loc$LOC, ": ", conditionMessage(e)); NULL
    })
    if (is.null(locus)) {
        for (j in seq_len(nrow(pairs)))
            results[[length(results) + 1]] <- result_row(loc, pairs$trait1[j], pairs$trait2[j],
                                                         "locus_not_processed")
        cat(sprintf("[%d/%d] locus %s: not processed by LAVA\n", n_done, n_todo, loc$LOC))
        next
    }
    univ <- run.univ(locus)
    univ_out[[length(univ_out) + 1]] <- data.frame(
        locus = as.character(loc$LOC), chr = loc$CHR, start = loc$START, stop = loc$STOP,
        phenotype = univ$phen, h2_obs = univ$h2.obs,
        h2_latent = if ("h2.latent" %in% names(univ)) univ$h2.latent else NA,
        p = univ$p, n_snps = locus$n.snps, K = locus$K,
        loaded_traits = paste(phenos, collapse = ";"), stringsAsFactors = FALSE)

    for (j in seq_len(nrow(pairs))) {
        t1 <- pairs$trait1[j]; t2 <- pairs$trait2[j]
        u1 <- univ[univ$phen == t1, , drop = FALSE]
        u2 <- univ[univ$phen == t2, , drop = FALSE]
        if (nrow(u1) == 0 || nrow(u2) == 0) {
            status <- "trait_dropped_in_locus"   # LAVA drops traits with invalid local estimates
            bv <- NULL
        } else if (is.na(u1$p) || is.na(u2$p) || u1$p >= univ_thresh || u2$p >= univ_thresh) {
            status <- "univ_filtered"
            bv <- NULL
        } else {
            set.seed(seed_for(base_seed, loc$LOC, t1, t2))
            bv <- tryCatch(
                run.bivar(locus, phenos = c(t1, t2), adap.thresh = adap_thresh,
                          p.values = do_p, CIs = do_ci, param.lim = param_lim),
                error = function(e) { message("locus ", loc$LOC, " ", t1, "-", t2, ": ",
                                              conditionMessage(e)); NULL })
            status <- if (is.null(bv)) "bivar_error"
                      else if (is.na(bv$rho[1])) "estimate_out_of_bounds"
                      else "ok"
        }
        results[[length(results) + 1]] <- result_row(loc, t1, t2, status,
                                                     n_snps = locus$n.snps, k = locus$K,
                                                     u1 = if (nrow(u1)) u1, u2 = if (nrow(u2)) u2,
                                                     bv = bv)
    }
    # progress line, visible in the task's .command.log while it runs
    cat(sprintf("[%d/%d] locus %s chr%s: %d SNPs, K=%d, %.1fs (elapsed %.1f min)\n",
                n_done, n_todo, loc$LOC, loc$CHR, locus$n.snps, locus$K,
                as.numeric(difftime(Sys.time(), t_locus, units = "secs")),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

res <- do.call(rbind, results)
write.table(res, paste0(out_prefix, ".lava.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
univ_df <- if (length(univ_out)) do.call(rbind, univ_out) else
    data.frame(locus = character(), chr = character(), start = numeric(), stop = numeric(),
               phenotype = character(), h2_obs = numeric(), h2_latent = numeric(),
               p = numeric(), n_snps = numeric(), K = numeric(), loaded_traits = character())
write.table(univ_df, paste0(out_prefix, ".lava_univ.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat("Status counts:\n"); print(table(res$status))
