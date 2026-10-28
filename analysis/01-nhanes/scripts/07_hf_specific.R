#!/usr/bin/env Rscript
# 07_hf_specific.R -- NHANES: is the bone-density association specific to heart
# failure rather than to cardiovascular disease in general?
#
# Motivation: the Mendelian randomisation layer found a BMD -> heart-failure effect
# and no BMD -> coronary effect. This script asks the same question in the population
# sample, with the same four adjustment sets used in 03_survey_analysis.R, and tests
# formally whether the BMD effect differs across the four cardiovascular outcomes.
#
# Design, weights and covariates are identical to 03_survey_analysis.R so the numbers
# are directly comparable: svydesign(ids=~psu, strata=~cycle:strata, weights=WTMEC2YR/5).
#
# Inputs : data/processed/nhanes/analysis_dataset.csv
# Outputs: results/tables/01-nhanes/
#            table10_hf_prevalence_by_bmd.csv
#            table11_hf_logistic.csv
#            table12_outcome_comparison.csv     head-to-head across outcomes
#            table13_outcome_heterogeneity.csv  Wald test for difference across outcomes
suppressPackageStartupMessages({
  library(survey)
})

tabdir <- "results/tables/01-nhanes"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv("data/processed/nhanes/analysis_dataset.csv", stringsAsFactors = FALSE)
d$bmd_cat <- factor(d$bmd_cat, levels = c("Normal", "Osteopenia", "Osteoporosis"))
d$sex <- factor(d$sex, levels = c("Male", "Female"))
d$race <- factor(d$race)
d$race <- relevel(d$race, ref = "NH White")
d$smoking <- factor(d$smoking, levels = c("Never", "Former", "Current"))
d$education <- factor(d$education, levels = c(">High school", "High school/GED", "<High school"))
d$cycle <- factor(d$cycle)
d$strata2 <- factor(paste(d$cycle, d$strata, sep = "_"))
d$wt5 <- d$wt_mec / 5
d$t_neg <- -d$t_score
d$bmd_low <- d$fn_bmd / 0.1

# four cardiovascular outcomes. `hf`, `chd`, `stroke` are already 0/1 numeric in the
# harmonised dataset; MI has to be built from the raw MCQ item, which the CSV stores as
# the strings "True"/"False".
d$chf <- as.numeric(d$hf)
d$mi <- as.numeric(d$mcq_mi == "True" | d$mcq_mi == TRUE)
d$cvd_death_num <- as.numeric(d$cvd_death)
d <- d[!is.na(d$chf), ]
stopifnot(sum(d$chf, na.rm = TRUE) > 100, sum(d$mi, na.rm = TRUE) > 100)

des <- svydesign(ids = ~psu, strata = ~strata2, weights = ~wt5, data = d, nest = TRUE)
outs <- c("chf", "chd", "mi", "stroke")

cov_std <- c("age", "sex", "race", "cycle")
cov_ses <- c("education", "pir", "bmi", "smoking")
cov_metab <- c("diabetes", "hypertension", "total_chol", "hdl", "egfr")
cov_bone <- c("albumin", "calcium", "phosphate", "alp")
models <- list(
  M1 = character(0),
  M2 = cov_std,
  M3 = c(cov_std, cov_ses, cov_metab),
  M4 = c(cov_std, cov_ses, cov_metab, cov_bone)
)

# ---------------------------------------------------------------- 1. prevalence
prev <- list()
for (o in outs) {
  f <- as.formula(paste("~", o))
  tab <- svyby(f, ~bmd_cat, des, svymean, na.rm = TRUE)
  # svyby names the estimate after the variable and the SE simply "se"
  est <- tab[[o]]
  sev <- tab[["se"]]
  prev[[o]] <- data.frame(outcome = o, bmd_cat = tab$bmd_cat,
                          prevalence = as.numeric(est), se = as.numeric(sev))
}
p1 <- do.call(rbind, prev)
# ratio osteoporosis vs normal, for a single comparable number per outcome
rat <- do.call(rbind, lapply(split(p1, p1$outcome), function(x) {
  get <- function(l) x$prevalence[x$bmd_cat == l]
  f <- get("Osteoporosis") / get("Normal")
  # delta-method SE on the log ratio
  lo <- x$se[x$bmd_cat == "Osteoporosis"] / get("Osteoporosis")
  ln <- x$se[x$bmd_cat == "Normal"] / get("Normal")
  data.frame(outcome = x$outcome[1], prev_normal = get("Normal"),
             prev_osteoporosis = get("Osteoporosis"), ratio = f,
             lo = exp(log(f) - 1.96 * sqrt(lo^2 + ln^2)),
             hi = exp(log(f) + 1.96 * sqrt(lo^2 + ln^2)))
}))
write.csv(p1, file.path(tabdir, "table10_hf_prevalence_by_bmd.csv"), row.names = FALSE)
cat("\nweighted prevalence by BMD category\n")
print(rat, row.names = FALSE)

# ---------------------------------------------------------------- 2. logistic
fit_or <- function(outcome, rhs, subset_expr = NULL) {
  vars <- unique(c(outcome, rhs, "bmd_cat", "t_neg", "bmd_low"))
  design <- if (is.null(subset_expr)) des else subset(des, subset_expr)
  res <- list()
  terms <- c("bmd_catOsteopenia", "bmd_catOsteoporosis")
  for (t in terms) {
    f <- as.formula(paste(outcome, "~", paste(c(rhs, "bmd_cat"), collapse = " + ")))
    fit <- try(svyglm(f, design = design, family = quasibinomial()), silent = TRUE)
    if (inherits(fit, "try-error")) next
    co <- summary(fit)$coefficients
    if (!t %in% rownames(co)) next
    b <- co[t, 1]; s <- co[t, 2]
    res[[t]] <- data.frame(outcome = outcome, term = t, or = exp(b),
                           lo = exp(b - 1.96 * s), hi = exp(b + 1.96 * s),
                           p = co[t, 4], stringsAsFactors = FALSE)
  }
  # per 1 SD lower BMD
  f <- as.formula(paste(outcome, "~", paste(c(rhs, "t_neg"), collapse = " + ")))
  fit <- try(svyglm(f, design = design, family = quasibinomial()), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    co <- summary(fit)$coefficients
    if ("t_neg" %in% rownames(co)) {
      b <- co["t_neg", 1]; s <- co["t_neg", 2]
      res[["t_neg"]] <- data.frame(outcome = outcome, term = "per SD lower T-score",
                                   or = exp(b), lo = exp(b - 1.96 * s),
                                   hi = exp(b + 1.96 * s), p = co["t_neg", 4],
                                   stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, res)
}

rows <- list()
for (o in outs) for (m in names(models)) {
  r <- fit_or(o, models[[m]])
  if (is.null(r) || !nrow(r)) next
  r$model <- m
  rows[[paste(o, m)]] <- r
}
p2 <- do.call(rbind, rows)
p2 <- p2[, c("outcome", "model", "term", "or", "lo", "hi", "p")]
write.csv(p2, file.path(tabdir, "table11_hf_logistic.csv"), row.names = FALSE)
cat("\nCHF vs other outcomes, fully adjusted (M4)\n")
print(p2[p2$model == "M4", ], row.names = FALSE)

# ---------------------------------------------------------------- 3. head to head
# one design per outcome, the same covariates, the same exposure term
head_to_head <- p2[p2$term == "per SD lower T-score", ]
p3 <- reshape(head_to_head[, c("outcome", "model", "or", "lo", "hi", "p")],
              idvar = "model", timevar = "outcome", direction = "wide")
write.csv(p3, file.path(tabdir, "table12_outcome_comparison.csv"), row.names = FALSE)
cat("\nper SD lower T-score, by outcome and adjustment set\n")
print(p3, row.names = FALSE)

# ---------------------------------------------------------------- 4. heterogeneity
# Wald test: is the T-score coefficient the same for CHF as for each other outcome?
# Fit the two outcomes jointly in one stacked design and test the interaction
# outcome x exposure. Stacking is done by duplicating rows with an outcome label and
# an indicator, which keeps the survey design valid (same PSUs, same weights/2).
het <- list()
for (o in setdiff(outs, "chf")) {
  long <- rbind(
    data.frame(d[, c("psu", "strata2", "wt5", "age", "sex", "race", "cycle", "education",
                     "pir", "bmi", "smoking", "diabetes", "hypertension", "total_chol",
                     "hdl", "egfr", "albumin", "calcium", "phosphate", "alp", "t_neg")],
               event = d$chf, which = "chf"),
    data.frame(d[, c("psu", "strata2", "wt5", "age", "sex", "race", "cycle", "education",
                     "pir", "bmi", "smoking", "diabetes", "hypertension", "total_chol",
                     "hdl", "egfr", "albumin", "calcium", "phosphate", "alp", "t_neg")],
               event = d[[o]], which = o))
  long$which <- factor(long$which, levels = c("chf", o))
  long <- long[!is.na(long$event), ]
  long$wt2 <- long$wt5 / 2
  # PSUs must be unique between the two copies or the design sees duplicated clusters
  long$psu2 <- paste(long$which, long$psu, sep = "_")
  des2 <- svydesign(ids = ~psu2, strata = ~strata2, weights = ~wt2, data = long, nest = TRUE)
  covs <- paste(c(cov_std, cov_ses, cov_metab, cov_bone), collapse = " + ")
  f <- as.formula(paste("event ~ t_neg * which +", covs))
  fit <- try(svyglm(f, design = des2, family = quasibinomial()), silent = TRUE)
  if (inherits(fit, "try-error")) next
  co <- summary(fit)$coefficients
  it <- grep("t_neg:which", rownames(co), value = TRUE)
  if (!length(it)) next
  z <- co[it[1], 1] / co[it[1], 2]
  het[[o]] <- data.frame(comparison = paste0("CHF vs ", o),
                         diff_in_per_SD = co[it[1], 1], se = co[it[1], 2],
                         z = z, p = 2 * pnorm(abs(z), lower.tail = FALSE))
}
if (length(het)) {
  p4 <- do.call(rbind, het)
  write.csv(p4, file.path(tabdir, "table13_outcome_heterogeneity.csv"), row.names = FALSE)
  cat("\nheterogeneity of the BMD effect: CHF versus each other outcome\n")
  print(p4, row.names = FALSE)
}

# ---------------------------------------------------------------- 5. prospective
# baseline BMD -> subsequent cardiovascular death, and CHF at baseline excluded
# so the association is temporally ordered
des_p <- subset(des, d$cvd == 0 | is.na(d$cvd))
rows <- list()
for (m in names(models)) {
  f <- as.formula(paste("cvd_death_num ~", paste(c(models[[m]], "t_neg"), collapse = " + ")))
  fit <- try(svyglm(f, design = des_p, family = quasibinomial()), silent = TRUE)
  if (inherits(fit, "try-error")) next
  co <- summary(fit)$coefficients
  rows[[m]] <- data.frame(model = m, n_events = sum(des_p$variables$cvd_death_num, na.rm = TRUE),
                          or = exp(co["t_neg", 1]),
                          lo = exp(co["t_neg", 1] - 1.96 * co["t_neg", 2]),
                          hi = exp(co["t_neg", 1] + 1.96 * co["t_neg", 2]),
                          p = co["t_neg", 4])
}
if (length(rows)) {
  p5 <- do.call(rbind, rows)
  write.csv(p5, file.path(tabdir, "table14_incident_cvd_death.csv"), row.names = FALSE)
  cat("\nprospective: baseline CVD excluded, BMD -> CVD death (per SD lower T-score)\n")
  print(p5, row.names = FALSE)
}

cat("\ndone\n")

