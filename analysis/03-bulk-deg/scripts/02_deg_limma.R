#!/usr/bin/env Rscript
# 02_deg_limma.R -- differential expression per GEO dataset with limma.
#
# Reads data/processed/geo/<GSE>/{expression.csv,phenotype.csv} produced by
# 01_series_matrix.R and writes results/tables/03-bulk-deg/deg_<GSE>.csv.
#
# Group assignment comes from analysis/03-bulk-deg/config/datasets.csv
# (columns: gse, case_regex, control_regex) matched against the sample title,
# source and characteristics fields. Datasets without a config row are skipped and
# their distinct characteristic values are printed so the config can be filled in
# from what is actually in the file, not from memory.
#
# Usage:
#   Rscript analysis/03-bulk-deg/scripts/02_deg_limma.R
#   Rscript analysis/03-bulk-deg/scripts/02_deg_limma.R GSE100927 GSE43292
#   Rscript analysis/03-bulk-deg/scripts/02_deg_limma.R --selftest

suppressPackageStartupMessages({
  library(limma)
  library(statmod)
})

script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(sub("^--file=", "", a[1])) else "analysis/03-bulk-deg/scripts"
}

read_dataset <- function(pr, gse) {
  dd <- file.path(pr, "data", "processed", "geo", gse)
  ef <- file.path(dd, "expression.csv")
  pf <- file.path(dd, "phenotype.csv")
  if (!file.exists(ef) || !file.exists(pf)) return(NULL)
  e <- as.matrix(read.csv(ef, row.names = 1, check.names = FALSE))
  p <- read.csv(pf, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(e) != nrow(p)) stop(gse, ": expression and phenotype sample counts differ")
  list(expr = e, pheno = p)
}

assign_groups <- function(p, case_regex, control_regex) {
  # match on the sample title and the characteristics block only: several series carry a
  # constant source_name such as "circulating monocytes in low-pbm phase", which would
  # otherwise assign every sample to the case group
  hay <- tolower(paste(p$title, p$characteristics, sep = " | "))
  hay[is.na(hay)] <- ""
  is_case <- grepl(tolower(case_regex), hay, perl = TRUE)
  is_ctrl <- grepl(tolower(control_regex), hay, perl = TRUE)
  g <- rep(NA_character_, length(hay))
  g[is_case & !is_ctrl] <- "case"
  g[is_ctrl & !is_case] <- "control"
  g
}

fit_limma <- function(e, group, block = NULL) {
  keep <- !is.na(group)
  e <- e[, keep, drop = FALSE]
  group <- factor(group[keep], levels = c("control", "case"))
  if (min(table(group)) < 2) stop("fewer than 2 samples in a group")
  design <- model.matrix(~group)
  if (!is.null(block)) {
    # paired samples (e.g. plaque and intact intima from the same patient): model the
    # donor as a random effect instead of treating the pairs as independent
    block <- factor(block[keep])
    corfit <- duplicateCorrelation(e, design, block = block)
    fit <- lmFit(e, design, block = block, correlation = corfit$consensus)
  } else {
    fit <- lmFit(e, design)
  }
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = "groupcase", number = Inf, sort.by = "P")
  data.frame(gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
             t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
             n_case = sum(group == "case"), n_control = sum(group == "control"),
             stringsAsFactors = FALSE)
}

pair_block <- function(p, field) {
  # characteristics look like "patient: 1 | disease state: hypertension | ..."
  if (!nzchar(field)) return(NULL)
  pat <- paste0("(?:^|\\| )", field, ": ([^|]+)")
  m <- regmatches(p$characteristics, regexpr(pat, p$characteristics, perl = TRUE))
  v <- trimws(sub(paste0("^(?:.*\\| )?", field, ": "), "", m))
  v[!nzchar(v)] <- NA
  if (all(is.na(v))) NULL else v
}

report_unconfigured <- function(gse, p) {
  vals <- unique(unlist(strsplit(paste(na.omit(p$characteristics), collapse = " | "), " \\| ")))
  message("  no config row: distinct characteristics -> ",
          paste(head(unique(trimws(vals)), 12), collapse = "; "))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  pr <- normalizePath(file.path(script_dir(), "..", "..", ".."))
  cfg_file <- file.path(pr, "analysis", "03-bulk-deg", "config", "datasets.csv")
  out_dir <- file.path(pr, "results", "tables", "03-bulk-deg")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  geo_dir <- file.path(pr, "data", "processed", "geo")
  gse <- args[!startsWith(args, "--")]
  if (!length(gse)) gse <- list.dirs(geo_dir, full.names = FALSE, recursive = FALSE)
  gse <- gse[!startsWith(gse, "_")]

  cfg <- if (file.exists(cfg_file)) read.csv(cfg_file, stringsAsFactors = FALSE) else
    data.frame(gse = character(0), case_regex = character(0), control_regex = character(0))
  rows <- list()
  for (g in gse) {
    d <- read_dataset(pr, g)
    if (is.null(d)) { message("skip ", g, " (not processed)"); next }
    row <- cfg[cfg$gse == g, ]
    if (!nrow(row)) { message("skip ", g); report_unconfigured(g, d$pheno); next }
    grp <- assign_groups(d$pheno, row$case_regex[1], row$control_regex[1])
    n_amb <- sum(grepl(tolower(row$case_regex[1]), tolower(paste(d$pheno$title, d$pheno$characteristics))) &
                 grepl(tolower(row$control_regex[1]), tolower(paste(d$pheno$title, d$pheno$characteristics))))
    message(sprintf("== %s: %d case / %d control / %d unassigned%s", g,
                    sum(grp == "case", na.rm = TRUE), sum(grp == "control", na.rm = TRUE),
                    sum(is.na(grp)), if (n_amb) sprintf(" (ambiguous: %d)", n_amb) else ""))
    if (sum(grp == "case", na.rm = TRUE) < 2 || sum(grp == "control", na.rm = TRUE) < 2) {
      message("  too few samples per group; skipped"); next
    }
    paired <- if ("paired_by" %in% names(row)) pair_block(d$pheno, as.character(row$paired_by[1])) else NULL
    if (!is.null(paired)) message(sprintf("  paired design on '%s' (%d levels)", row$paired_by[1],
                                          length(unique(na.omit(paired[!is.na(grp)])))))
    tt <- fit_limma(d$expr, grp, if (is.null(paired)) NULL else paired)
    f <- file.path(out_dir, paste0("deg_", g, ".csv"))
    write.csv(tt, f, row.names = FALSE)
    n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 1)
    n_dn <- sum(tt$adj.P.Val < 0.05 & tt$logFC < -1)
    n_sig <- sum(tt$adj.P.Val < 0.05)
    message(sprintf("  -> %s: %d genes, FDR<0.05 %d (%d up / %d down at |logFC|>1)",
                    basename(f), nrow(tt), n_sig, n_up, n_dn))
    rows[[g]] <- data.frame(gse = g, genes = nrow(tt), fdr05 = n_sig, up = n_up, down = n_dn,
                            n_case = sum(grp == "case", na.rm = TRUE),
                            n_control = sum(grp == "control", na.rm = TRUE),
                            paired = !is.null(paired))
  }
  if (length(rows)) {
    write.csv(do.call(rbind, rows), file.path(out_dir, "deg_summary.csv"), row.names = FALSE)
  } else {
    message("nothing analysed (no processed datasets or no config rows)")
  }
}

selftest <- function() {
  set.seed(42)
  ng <- 800
  n <- 12
  e <- matrix(rnorm(ng * n, 8, 1), ng, dimnames = list(sprintf("G%04d", seq_len(ng)),
                                                       sprintf("S%02d", seq_len(n))))
  e[1:40, 7:12] <- e[1:40, 7:12] + 2          # 40 true up-regulated genes
  group <- rep(c("control", "case"), each = 6)
  tt <- fit_limma(e, group)
  true_up <- sprintf("G%04d", 1:40)
  rank_of <- match(true_up, tt$gene)
  in_top10pct <- sum(rank_of <= ng / 10)
  cat(sprintf("selftest: %d/40 planted genes within the top 10%% of the ranking, median logFC %.2f\n",
              in_top10pct, median(tt$logFC[1:40])))
  pass <- in_top10pct >= 35 && median(tt$logFC[1:40]) > 1.5
  cat(if (pass) "PASS\n" else "FAIL\n")
  if (!pass) quit(status = 1)
}

if ("--selftest" %in% commandArgs(trailingOnly = TRUE)) selftest() else main()

