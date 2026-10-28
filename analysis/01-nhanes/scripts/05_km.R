# Step 5 - survey-weighted Kaplan-Meier survival by femoral neck BMD category.
#
# survey::svykm() with design-based standard errors did not finish on this sample
# (n = 14,260, 5 pooled cycles), so the curves are computed directly from the MEC
# weights (Horvitz-Thompson risk sets) with weighted Greenwood variance. They are
# descriptive; inference for the BMD-mortality association comes from the
# design-based Cox models in 03_survey_analysis.R (Table 4).
#
# Input : data/processed/nhanes/analysis_dataset.csv
# Output: results/tables/01-nhanes/table9_km_allcause_by_bmd.csv
suppressPackageStartupMessages(library(survival))

tabdir <- "results/tables/01-nhanes"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv("data/processed/nhanes/analysis_dataset.csv", stringsAsFactors = FALSE)
d$bmd_cat <- factor(d$bmd_cat, levels = c("Normal", "Osteopenia", "Osteoporosis"))
d <- d[!is.na(d$fu_years_exm) & !is.na(d$MORTSTAT) & !is.na(d$bmd_cat) & d$fu_years_exm > 0, ]
d$wt5 <- d$wt_mec / 5
# Normalise weights to mean 1 so the weighted risk sets stay on the participant
# scale and the Greenwood variance is not collapsed by the raw survey weights.
d$wn <- d$wt5 / mean(d$wt5, na.rm = TRUE)
d$event <- as.numeric(d$MORTSTAT)

weighted_km <- function(time, event, weight) {
  ord <- order(time)
  time <- time[ord]; event <- event[ord]; weight <- weight[ord]
  ut <- sort(unique(time[event == 1]))
  n_risk <- sapply(ut, function(t) sum(weight[time >= t]))
  d_ev <- sapply(ut, function(t) sum(weight[time == t & event == 1]))
  surv <- cumprod(1 - d_ev / n_risk)
  # weighted Greenwood: var(S) = S^2 * sum(d / (n(n-d)))
  gr <- cumsum(d_ev / (n_risk * pmax(n_risk - d_ev, 1e-12)))
  se <- surv * sqrt(gr)
  data.frame(time = c(0, ut), surv = c(1, surv), se = c(0, se))
}

grid_t <- seq(0, 14, by = 0.25)
res <- do.call(rbind, lapply(levels(d$bmd_cat), function(g) {
  sub <- d[d$bmd_cat == g, ]
  k <- weighted_km(sub$fu_years_exm, sub$event, sub$wn)
  idx <- findInterval(grid_t, k$time, left.open = FALSE)
  s <- k$surv[pmax(idx, 1)]
  se <- k$se[pmax(idx, 1)]
  data.frame(stratum = g, time = grid_t, surv = s,
             lcl = pmax(0, s - 1.96 * se), ucl = pmin(1, s + 1.96 * se),
             n_unweighted = nrow(sub),
             deaths_unweighted = sum(sub$event))
}))

write.csv(res, file.path(tabdir, "table9_km_allcause_by_bmd.csv"), row.names = FALSE)
cat("strata:", levels(d$bmd_cat), "\n")
print(res[res$time %in% c(0, 5, 10, 14), c("stratum", "time", "surv", "lcl", "ucl")],
      row.names = FALSE, digits = 3)
cat("->", file.path(tabdir, "table9_km_allcause_by_bmd.csv"), "\n")

