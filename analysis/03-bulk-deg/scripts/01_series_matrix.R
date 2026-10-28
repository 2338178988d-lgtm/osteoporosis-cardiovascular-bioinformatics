#!/usr/bin/env Rscript
# 01_series_matrix.R -- read GEO series matrix files, QC, normalise, collapse probes.
#
# Usage:
#   Rscript analysis/03-bulk-deg/scripts/01_series_matrix.R                 # all present
#   Rscript analysis/03-bulk-deg/scripts/01_series_matrix.R GSE100927 GSE43292
#   Rscript analysis/03-bulk-deg/scripts/01_series_matrix.R --selftest
#
# Output per dataset in data/processed/geo/<GSE>/: expression.csv (genes x samples,
# log2), phenotype.csv, qc.csv.

suppressPackageStartupMessages(library(limma))

parse_series_matrix <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  meta <- list()
  rows <- character(0)
  in_table <- FALSE
  repeat {
    ln <- readLines(con, n = 1, warn = FALSE)
    if (!length(ln)) break
    if (grepl("^!series_matrix_table_begin", ln)) { in_table <- TRUE; next }
    if (grepl("^!series_matrix_table_end", ln)) break
    if (in_table) {
      rows <- c(rows, ln)
    } else if (grepl("^!", ln)) {
      key <- sub("^!", "", strsplit(ln, "\t", fixed = TRUE)[[1]][1])
      vals <- gsub('^"|"$', "", strsplit(ln, "\t", fixed = TRUE)[[1]][-1])
      meta[[key]] <- c(meta[[key]], list(vals))
    }
  }
  if (!length(rows)) stop("no series matrix table in ", path)
  tab <- read.delim(text = paste(rows, collapse = "\n"), check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "")
  list(meta = meta, table = tab)
}

sample_meta <- function(meta, tab) {
  n <- ncol(tab) - 1
  pick <- function(k) {
    v <- meta[[k]]
    if (is.null(v)) return(rep(NA_character_, n))
    v <- v[[length(v)]]
    if (length(v) == n) v else rep(NA_character_, n)
  }
  chars <- rep(NA_character_, n)
  for (v in meta[["Sample_characteristics_ch1"]]) {
    if (length(v) == n) chars <- ifelse(is.na(chars), v, paste(chars, v, sep = " | "))
  }
  data.frame(sample = pick("Sample_geo_accession"),
             title = pick("Sample_title"),
             source = pick("Sample_source_name_ch1"),
             characteristics = chars,
             platform = pick("Sample_platform_id"),
             stringsAsFactors = FALSE)
}

assay_numeric <- function(tab) {
  m <- as.matrix(tab[, -1, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "double")
  # the matrix header and ID column are quoted in the file; quote = "" keeps the
  # quote characters, so strip them from the identifiers
  rownames(m) <- gsub('^"|"$', "", as.character(tab[[1]]))
  colnames(m) <- gsub('^"|"$', "", colnames(m))
  m[is.finite(rowSums(m)), , drop = FALSE]
}

looks_logged <- function(m) {
  q <- quantile(m, c(0.01, 0.99), na.rm = TRUE)
  q[2] < 30 && q[1] >= 0
}

collapse_to_gene <- function(m, symbols) {
  if (is.null(symbols) || all(is.na(symbols))) return(list(expr = m, level = "probe"))
  keep <- !is.na(symbols) & symbols != ""
  m <- m[keep, , drop = FALSE]
  symbols <- symbols[keep]
  ord <- order(rowMeans(m), decreasing = TRUE)
  m <- m[ord, , drop = FALSE]
  symbols <- symbols[ord]
  dup <- duplicated(symbols)
  expr <- m[!dup, , drop = FALSE]
  rownames(expr) <- symbols[!dup]
  list(expr = expr, level = "gene")
}

read_symbols <- function(gse, platform_dir) {
  f <- file.path(platform_dir, paste0(gse, "_symbols.tsv"))
  if (!file.exists(f)) return(NULL)
  d <- read.delim(f, stringsAsFactors = FALSE)
  setNames(d$symbol, d$probe)
}

run_dataset <- function(gse, raw_dir, out_dir, platform_dir) {
  f <- file.path(raw_dir, gse, paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(f)) f <- file.path(raw_dir, gse, paste0(gse, "_series_matrix.txt"))
  if (!file.exists(f)) { message("skip ", gse, " (not downloaded)"); return(FALSE) }
  message("== ", gse)
  parsed <- parse_series_matrix(f)
  tab <- parsed$table
  m <- tryCatch(assay_numeric(tab), error = function(e) { message("  ", conditionMessage(e)); NULL })
  if (is.null(m) || !ncol(m)) { message("  no numeric matrix"); return(FALSE) }
  ph <- sample_meta(parsed$meta, tab)
  logged <- looks_logged(m)
  message(sprintf("  %d probes x %d samples, values %s", nrow(m), ncol(m),
                  if (logged) "already log2" else "linear -> log2"))
  if (!logged) m <- log2(pmax(m, 1))
  m <- normalizeBetweenArrays(m, method = "quantile")
  sym <- read_symbols(gse, platform_dir)
  mapped <- if (is.null(sym)) 0 else sum(!is.na(sym[rownames(m)]))
  message(sprintf("  probes with a gene symbol: %d of %d", mapped, nrow(m)))
  res <- collapse_to_gene(m, if (is.null(sym)) NULL else unname(sym[rownames(m)]))
  od <- file.path(out_dir, gse)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  write.csv(round(res$expr, 5), file.path(od, "expression.csv"), quote = FALSE)
  write.csv(ph, file.path(od, "phenotype.csv"), row.names = FALSE)
  write.csv(data.frame(sample = ph$sample,
                       median = apply(m, 2, median),
                       iqr = apply(m, 2, IQR),
                       pct_zero = apply(m, 2, function(x) mean(x == 0))),
            file.path(od, "qc.csv"), row.names = FALSE)
  message(sprintf("  -> %s: %d %ss x %d samples", od, nrow(res$expr), res$level, ncol(res$expr)))
  TRUE
}

selftest <- function() {
  tmp <- tempfile("geo_selftest_")
  dir.create(file.path(tmp, "raw", "GSE000001"), recursive = TRUE)
  dir.create(file.path(tmp, "out"), recursive = TRUE)
  set.seed(1)
  nprobe <- 500
  nsamp <- 8
  base <- 2^rnorm(nprobe, 8, 1.5)
  lin <- sapply(seq_len(nsamp), function(i) base * exp(rnorm(nprobe, 0, 0.25)))
  lin[1:100, 5:8] <- lin[1:100, 5:8] * 3      # planted 3-fold effect in 100 probes
  tab <- data.frame(ID_REF = sprintf("ILMN_%05d", seq_len(nprobe)), lin, check.names = FALSE)
  head <- c(
    "!Series_title\t\"synthetic\"",
    paste0("!Sample_geo_accession\t", paste(sprintf('"GSM%04d"', 1:nsamp), collapse = "\t")),
    paste0("!Sample_title\t", paste(sprintf('"s%d"', 1:nsamp), collapse = "\t")),
    paste0("!Sample_characteristics_ch1\t", paste(c(rep('"group: ctrl"', 4), rep('"group: case"', 4)),
                                                  collapse = "\t")),
    "!series_matrix_table_begin")
  f <- file.path(tmp, "raw", "GSE000001", "GSE000001_series_matrix.txt.gz")
  con <- gzfile(f, "w")
  writeLines(c(head, paste(names(tab), collapse = "\t")), con)
  write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines("!series_matrix_table_end", con)
  close(con)

  ok <- run_dataset("GSE000001", file.path(tmp, "raw"), file.path(tmp, "out"), tmp)
  e <- as.matrix(read.csv(file.path(tmp, "out", "GSE000001", "expression.csv"),
                          row.names = 1, check.names = FALSE))
  target <- as.matrix(read.csv(file.path(tmp, "out", "GSE000001", "phenotype.csv")))
  med <- apply(e, 2, median)
  spread <- max(med) - min(med)
  cat(sprintf("\nselftest: %d x %d, median range across samples %.3f (quantile normalised)\n",
              nrow(e), ncol(e), spread))
  cat(sprintf("selftest: all values finite = %s\n", all(is.finite(e))))
  pass <- ok && spread < 0.01 && all(is.finite(e)) && ncol(e) == nsamp
  cat(if (pass) "PASS\n" else "FAIL\n")
  unlink(tmp, recursive = TRUE)
  if (!pass) quit(status = 1)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--selftest" %in% args) return(selftest())
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1])) else "analysis/03-bulk-deg/scripts"
  pr <- normalizePath(file.path(here, "..", "..", ".."))
  raw <- file.path(pr, "data", "raw", "geo")
  out <- file.path(pr, "data", "processed", "geo")
  plat <- file.path(pr, "data", "processed", "geo", "_platform")
  dir.create(plat, recursive = TRUE, showWarnings = FALSE)
  gse <- args[!startsWith(args, "--")]
  if (!length(gse)) {
    gse <- sub("_series_matrix\\.txt(\\.gz)?$", "",
               list.files(raw, recursive = TRUE, pattern = "_series_matrix\\.txt(\\.gz)?$"))
    gse <- unique(basename(gse))
  }
  if (!length(gse)) { message("no series matrices under ", raw); return(invisible(NULL)) }
  ok <- vapply(gse, run_dataset, logical(1), raw_dir = raw, out_dir = out, platform_dir = plat)
  message(sprintf("\nprocessed %d/%d datasets", sum(ok), length(ok)))
}

main()

