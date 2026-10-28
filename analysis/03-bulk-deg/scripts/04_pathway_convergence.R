#!/usr/bin/env Rscript
# 04_pathway_convergence.R -- gene-level intersection of bone and vascular DEGs fails
# (03_intersect.R shows the overlap is at or below the chance expectation), so the
# question is asked one level up: do the two sides converge on the same *pathways*?
#
# For each dataset: ORA (clusterProfiler) against MSigDB Hallmark, restricted to the
# genes measurable on both tissue sides, so the two sides are compared on a common
# universe. Then the bone-side and vascular-side enriched sets are intersected and the
# observed count is compared with the count expected from random gene sets of the same
# size (permutation), which is what makes "convergence" a testable claim.
#
# Usage:
#   Rscript analysis/03-bulk-deg/scripts/04_pathway_convergence.R
#   Rscript analysis/03-bulk-deg/scripts/04_pathway_convergence.R --selftest
#
# Outputs (results/tables/03-bulk-deg/):
#   pathway_ora_<dataset>.csv     ORA result per dataset
#   pathway_convergence.csv       bone-vs-vascular shared pathways with permutation p
#   pathway_summary.csv           enriched pathway counts per dataset

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(stats)
})

script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(sub("^--file=", "", a[1])) else "analysis/03-bulk-deg/scripts"
}

hallmark <- function() {
  h <- msigdbr(species = "Homo sapiens", collection = "H")
  split(h$gene_symbol, h$gs_name)
}

read_deg <- function(dir, gse) {
  f <- file.path(dir, paste0("deg_", gse, ".csv"))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}

# ORA on the genes that are measurable on both sides, so the background is shared
ora_one <- function(d, universe, sets, fdr = 0.05, lfc = 0.0, direction = "both") {
  sig <- d[d$adj.P.Val < fdr, ]
  if (lfc > 0) sig <- sig[abs(sig$logFC) > lfc, ]
  if (direction == "up") sig <- sig[sig$logFC > 0, ]
  if (direction == "down") sig <- sig[sig$logFC < 0, ]
  genes <- intersect(sig$gene, universe)
  if (length(genes) < 5) return(NULL)
  sets_u <- lapply(sets, function(s) intersect(s, universe))
  sets_u <- sets_u[lengths(sets_u) >= 10 & lengths(sets_u) <= 500]
  if (!length(sets_u)) return(NULL)
  res <- try(enricher(genes, TERM2GENE = data.frame(
    term = rep(names(sets_u), lengths(sets_u)),
    gene = unlist(sets_u, use.names = FALSE)), universe = universe,
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res) || !nrow(as.data.frame(res))) return(NULL)
  out <- as.data.frame(res)
  out$n_sig_genes <- length(genes)
  out[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "n_sig_genes")]
}

perm_shared <- function(bone_sets, vasc_sets, all_sets, n_perm = 2000, seed = 1) {
  observed <- length(intersect(bone_sets, vasc_sets))
  names_all <- names(all_sets)
  sizes <- c(length(bone_sets), length(vasc_sets))
  # sampling without replacement needs each set to be selectable from the pool
  if (any(sizes > length(names_all))) {
    return(c(observed = observed, expected = NA, p = NA))
  }
  set.seed(seed)
  null <- replicate(n_perm, {
    a <- sample(names_all, sizes[1])
    b <- sample(names_all, sizes[2])
    length(intersect(a, b))
  })
  c(observed = observed, expected = mean(null), p = (sum(null >= observed) + 1) / (n_perm + 1))
}

main <- function() {
  pr <- normalizePath(file.path(script_dir(), "..", "..", ".."))
  tab <- file.path(pr, "results", "tables", "03-bulk-deg")
  bone <- c("GSE56815", "GSE56814", "GSE7158", "GSE7429", "GSE13850", "GSE35956")
  vasc <- c("GSE100927", "GSE43292", "GSE28829")
  sets <- hallmark()
  message("Hallmark gene sets: ", length(sets))

  degs <- list()
  for (g in c(bone, vasc, "GSE37558")) {
    d <- read_deg(tab, g)
    if (!is.null(d)) degs[[g]] <- d
  }
  if (!length(degs)) { message("no DEG tables found"); return(invisible(NULL)) }
  universe <- Reduce(intersect, lapply(degs, function(d) d$gene))
  message("shared measured universe: ", length(universe), " genes")

  ora <- list()
  summ <- list()
  for (g in names(degs)) {
    for (dir in c("up", "down")) {
      r <- ora_one(degs[[g]], universe, sets, direction = dir)
      if (is.null(r)) next
      r$dataset <- g
      r$direction <- dir
      r$side <- if (g %in% bone) "bone" else if (g %in% vasc) "vascular" else "mechanism"
      ora[[paste(g, dir)]] <- r
      sig <- r[r$p.adjust < 0.05, ]
      if (nrow(sig)) {
        f <- file.path(tab, paste0("pathway_ora_", g, "_", dir, ".csv"))
        write.csv(sig, f, row.names = FALSE)
      }
      summ[[paste(g, dir)]] <- data.frame(dataset = g, side = r$side[1], direction = dir,
                                          n_sig_genes = r$n_sig_genes[1],
                                          n_pathways_fdr05 = nrow(sig),
                                          min_padj = if (nrow(sig)) min(sig$p.adjust) else NA)
    }
  }
  if (length(summ)) {
    s <- do.call(rbind, summ)
    write.csv(s, file.path(tab, "pathway_summary.csv"), row.names = FALSE)
    message("\nenriched Hallmark pathways per dataset/direction:")
    print(s, row.names = FALSE)
  }

  # convergence: shared enriched pathways between the two sides
  side_sets <- function(side) {
    d <- unique(unlist(lapply(ora, function(r) {
      if (r$side[1] != side) return(character(0))
      r$ID[r$p.adjust < 0.05]
    })))
    d
  }
  bone_set <- side_sets("bone")
  vasc_set <- side_sets("vascular")
  mech_set <- side_sets("mechanism")
  ps <- perm_shared(bone_set, vasc_set, sets)
  conv <- data.frame(
    comparison = c("bone vs vascular", "bone vs mechanism", "vascular vs mechanism"),
    bone_like = c(length(bone_set), length(bone_set), length(vasc_set)),
    vascular_like = c(length(vasc_set), length(mech_set), length(mech_set)),
    observed = c(ps["observed"], perm_shared(bone_set, mech_set, sets)["observed"],
                 perm_shared(vasc_set, mech_set, sets)["observed"]),
    expected = c(ps["expected"], perm_shared(bone_set, mech_set, sets)["expected"],
                 perm_shared(vasc_set, mech_set, sets)["expected"]),
    perm_p = c(ps["p"], perm_shared(bone_set, mech_set, sets)["p"],
               perm_shared(vasc_set, mech_set, sets)["p"]))
  write.csv(conv, file.path(tab, "pathway_convergence.csv"), row.names = FALSE)
  message("\npathway convergence (permutation over Hallmark set names):")
  print(conv, row.names = FALSE)
  message("\nshared bone-vascular pathways: ",
          if (length(intersect(bone_set, vasc_set))) paste(intersect(bone_set, vasc_set), collapse = ", ")
          else "(none)")
}

selftest <- function() {
  # a synthetic universe with planted pathways; ORA must recover the planted one and the
  # permutation test must call a genuine overlap significant and a disjoint one not
  set.seed(9)
  universe <- sprintf("G%04d", 1:800)
  planted <- universe[1:60]
  sets <- list(PATH_A = planted, PATH_B = universe[101:160], PATH_C = universe[201:260],
               PATH_D = universe[301:360], PATH_E = universe[401:460],
               PATH_F = universe[501:560], PATH_G = universe[601:660])
  d <- data.frame(gene = universe,
                  logFC = ifelse(universe %in% planted, 1.6, rnorm(800, 0, 0.2)),
                  P.Value = ifelse(universe %in% planted, 1e-6, runif(800)),
                  adj.P.Val = ifelse(universe %in% planted, 1e-6, runif(800)),
                  stringsAsFactors = FALSE)
  r <- ora_one(d, universe, sets, direction = "up")
  ok1 <- !is.null(r) && r$ID[1] == "PATH_A"
  # realistic sizes: what the real analysis looks like is a few dozen enriched pathways
  # out of the 50 Hallmark sets, so the self-test uses that scale rather than 2-of-7,
  # where the permutation test has almost no power by construction
  big <- setNames(lapply(1:50, function(i) universe[((i * 10) %% 700 + 1):(((i * 10) %% 700) + 41)]),
                  sprintf("SET_%02d", 1:50))
  bone_hits <- sprintf("SET_%02d", c(1:10))
  vasc_hits <- sprintf("SET_%02d", c(1:8, 21:24))
  p_true <- perm_shared(bone_hits, vasc_hits, big, n_perm = 1000)
  p_null <- perm_shared(sprintf("SET_%02d", 31:40), sprintf("SET_%02d", 41:50), big, n_perm = 1000)
  cat(sprintf("selftest: top ORA pathway = %s (planted PATH_A)\n", if (is.null(r)) "none" else r$ID[1]))
  cat(sprintf("selftest: overlap of 10 vs 12 sets sharing 8: observed=%d expected=%.2f p=%.4f\n",
              p_true["observed"], p_true["expected"], p_true["p"]))
  cat(sprintf("selftest: disjoint 10 vs 10 sets: observed=%d expected=%.2f p=%.3f\n",
              p_null["observed"], p_null["expected"], p_null["p"]))
  pass <- ok1 && p_true["p"] < 0.05 && p_null["p"] > 0.2
  cat(if (pass) "PASS\n" else "FAIL\n")
  if (!pass) quit(status = 1)
}

if ("--selftest" %in% commandArgs(trailingOnly = TRUE)) selftest() else main()

