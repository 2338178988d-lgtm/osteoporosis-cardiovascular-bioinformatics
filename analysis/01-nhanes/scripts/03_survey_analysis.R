# Step 3 - survey-weighted analysis of the femoral neck BMD x CVD association
# and of BMD with subsequent all-cause / cardiovascular mortality.
#
# Design: NHANES 2005-2006, 2007-2008, 2009-2010, 2013-2014, 2017-2018 pooled.
# Weights: WTMEC2YR / 5 (five 2-year cycles). Strata nested within cycle.
# Subgroups use domain analysis (subset on the design), never pre-filtered data.
#
# Inputs : data/processed/nhanes/analysis_dataset.csv
# Outputs: results/tables/01-nhanes/*.csv
suppressPackageStartupMessages({
  library(survey)
  library(survival)
})

tabdir <- "results/tables/01-nhanes"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv("data/processed/nhanes/analysis_dataset.csv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------- preparation
d$bmd_cat <- factor(d$bmd_cat, levels = c("Normal", "Osteopenia", "Osteoporosis"))
d$bmd_ord <- as.numeric(d$bmd_cat)
d$sex <- factor(d$sex, levels = c("Male", "Female"))
d$race <- factor(d$race)
d$race <- relevel(d$race, ref = "NH White")
d$smoking <- factor(d$smoking, levels = c("Never", "Former", "Current"))
d$education <- factor(d$education, levels = c(">High school", "High school/GED", "<High school"))
d$cycle <- factor(d$cycle)
d$strata2 <- factor(paste(d$cycle, d$strata, sep = "_"))
d$wt5 <- d$wt_mec / 5
d$age65 <- factor(ifelse(d$age >= 65, "65+", "40-64"), levels = c("40-64", "65+"))
# OR > 1 means lower BMD -> higher risk
d$t_neg <- -d$t_score
d$bmd_low <- d$fn_bmd / 0.1

d$cvd <- as.numeric(d$cvd)
d$chd <- as.numeric(d$chd)
d$stroke <- as.numeric(d$stroke)
d$self_op <- as.numeric(d$self_op)
d$cvd_death_num <- as.numeric(d$cvd_death)
d$dead <- as.numeric(d$MORTSTAT)
stopifnot(sum(!is.na(d$cvd)) > 0, sum(!is.na(d$cvd_death_num)) > 0)

des <- svydesign(ids = ~psu, strata = ~strata2, weights = ~wt5, data = d, nest = TRUE)

cov_std <- c("age", "race", "education", "pir", "bmi", "smoking")
cov_metab <- c("diabetes", "hypertension", "total_chol", "hdl", "egfr")
cov_bone <- c("albumin", "calcium", "phosphate", "alp")
cov_all <- c(cov_std, cov_metab, cov_bone)

f <- function(rhs) as.formula(paste("~", paste(rhs, collapse = " + ")))

# ---------------------------------------------------------------- Table 1
# Weighted characteristics by femoral neck BMD category (WHO T-score).
t1 <- data.frame()
for (v in c("age", "fn_bmd", "t_score", "bmi", "pir")) {
  est <- svyby(as.formula(paste0("~", v)), ~bmd_cat, des, svymean, na.rm = TRUE, keep.names = FALSE)
  ci <- confint(est)
  t1 <- rbind(t1, data.frame(variable = v, level = "mean", bmd_cat = est$bmd_cat,
                             weighted_estimate = est[[v]], lcl = ci[, 1], ucl = ci[, 2],
                             n_unweighted = as.numeric(table(d$bmd_cat))))
}
cat_vars <- c("sex", "cvd", "chd", "stroke", "hypertension", "diabetes", "smoking",
              "education", "self_op", "race")
for (v in cat_vars) {
  est <- svyby(as.formula(paste0("~", v)), ~bmd_cat, des, svymean, na.rm = TRUE, keep.names = FALSE)
  for (nm in setdiff(names(est), c("bmd_cat", "se"))) {
    t1 <- rbind(t1, data.frame(variable = v, level = nm, bmd_cat = est$bmd_cat,
                               weighted_estimate = est[[nm]], lcl = NA, ucl = NA,
                               n_unweighted = NA))
  }
}
write.csv(t1, file.path(tabdir, "table1_weighted_by_bmd_category.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Table 2
prev <- svyby(~cvd, ~bmd_cat + sex, des, svymean, na.rm = TRUE, keep.names = FALSE, vartype = "ci")
prev_all <- svyby(~cvd, ~bmd_cat, des, svymean, na.rm = TRUE, keep.names = FALSE, vartype = "ci")
prev_all$sex <- "Both"
cols <- intersect(c("bmd_cat", "sex", "cvd", "ci_l", "ci_u"), names(prev_all))
prev_out <- rbind(prev_all[, cols], prev[, cols])
n_by <- as.data.frame(table(d$bmd_cat, d$sex))
names(n_by) <- c("bmd_cat", "sex", "n_unweighted")
prev_out <- merge(prev_out, n_by, by = c("bmd_cat", "sex"), all.x = TRUE)
names(prev_out)[names(prev_out) == "cvd"] <- "weighted_prevalence"
write.csv(prev_out, file.path(tabdir, "table2_cvd_prevalence_by_bmd.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Table 3: logistic
add <- function(a, b) if (is.null(a)) b else rbind(a, b)

extract_glm <- function(model, label, subgroup, exposure, outcome) {
  s <- summary(model)$coefficients
  ci <- suppressMessages(confint(model))
  keep <- grep("^bmd_cat|^t_neg|^bmd_low|^bmd_ord", rownames(s))
  if (length(keep) == 0) return(NULL)
  data.frame(outcome = outcome, model = label, subgroup = subgroup, exposure = exposure,
             term = rownames(s)[keep],
             OR = exp(s[keep, 1]), lcl = exp(ci[keep, 1]), ucl = exp(ci[keep, 2]),
             p = s[keep, 4], n = nrow(model$model), row.names = NULL)
}

run_logit <- function(outcome, design, subgroup = "All", drop = character(0)) {
  res <- NULL
  rhs_models <- list(
    "M1 crude" = c("bmd_cat"),
    "M2 +age/sex/race/cycle" = c("bmd_cat", "age", "sex", "race", "cycle"),
    "M3 +lifestyle/metabolic" = c("bmd_cat", "age", "sex", "race", "cycle", cov_std, cov_metab),
    "M4 +mineral metabolism" = c("bmd_cat", "age", "sex", "race", "cycle", cov_all)
  )
  rhs_models <- lapply(rhs_models, function(v) setdiff(v, drop))
  for (nm in names(rhs_models)) {
    form <- as.formula(paste(outcome, "~", paste(rhs_models[[nm]], collapse = " + ")))
    fit <- try(svyglm(form, design = design, family = quasibinomial()), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      res <- add(res, extract_glm(fit, nm, subgroup, "bmd_cat", outcome))
    }
  }
  # continuous exposures, fully adjusted
  for (expo in c("t_neg", "bmd_low")) {
    form <- as.formula(paste(outcome, "~", expo, "+",
                             paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + ")))
    fit <- try(svyglm(form, design = design, family = quasibinomial()), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      res <- add(res, extract_glm(fit, "M4 continuous", subgroup, expo, outcome))
    }
  }
  # p for trend across WHO categories (fully adjusted)
  form <- as.formula(paste(outcome, "~ bmd_ord +",
                           paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + ")))
  fit <- try(svyglm(form, design = design, family = quasibinomial()), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)$coefficients
    res <- add(res, data.frame(outcome = outcome, model = "M4 p-for-trend", subgroup = subgroup,
                               exposure = "bmd_ord", term = "bmd_ord",
                               OR = exp(s["bmd_ord", 1]), lcl = NA, ucl = NA,
                               p = s["bmd_ord", 4], n = nrow(fit$model)))
  }
  res
}

t3 <- rbind(
  run_logit("cvd", des),
  run_logit("cvd", subset(des, sex == "Male"), "Male", drop = "sex"),
  run_logit("cvd", subset(des, sex == "Female"), "Female", drop = "sex"),
  run_logit("cvd", subset(des, age65 == "40-64"), "Age 40-64", drop = "age"),
  run_logit("cvd", subset(des, age65 == "65+"), "Age 65+", drop = "age"),
  run_logit("chd", des),
  run_logit("stroke", des)
)
write.csv(t3, file.path(tabdir, "table3_survey_logistic_cvd.csv"), row.names = FALSE)
# ---------------------------------------------------------------- Table 4: Cox
run_cox <- function(time, event, design, subgroup = "All", drop = character(0)) {
  res <- list()
  base <- c("age", "sex", "race", "cycle")
  rhs_models <- list(
    "M1 crude" = c("bmd_cat"),
    "M2 +age/sex/race/cycle" = c("bmd_cat", base),
    "M3 +lifestyle/metabolic/prevalent CVD" = c("bmd_cat", base, cov_std, cov_metab, "cvd")
  )
  rhs_models <- lapply(rhs_models, function(v) setdiff(v, drop))
  for (nm in names(rhs_models)) {
    form <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ",
                              paste(rhs_models[[nm]], collapse = " + ")))
    fit <- try(svycoxph(form, design = design), silent = TRUE)
    if (inherits(fit, "try-error")) next
    s <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))
    keep <- grep("^bmd_cat", rownames(s))
    if (length(keep) == 0) next
    res[[length(res) + 1]] <- data.frame(outcome = event, model = nm, subgroup = subgroup,
                                         exposure = "bmd_cat", term = rownames(s)[keep],
                                         HR = s[keep, "exp(coef)"],
                                         lcl = exp(ci[keep, 1]), ucl = exp(ci[keep, 2]),
                                         p = s[keep, "Pr(>|z|)"], n = fit$n)
  }
  do.call(rbind, res)
}

t4 <- rbind(
  run_cox("fu_years_exm", "dead", des),
  run_cox("fu_years_exm", "cvd_death_num", des),
  run_cox("fu_years_exm", "dead", subset(des, sex == "Female"), "Female", drop = "sex"),
  run_cox("fu_years_exm", "dead", subset(des, sex == "Male"), "Male", drop = "sex")
)
write.csv(t4, file.path(tabdir, "table4_survey_cox_mortality.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Table 5: sensitivity
sens <- list()
# (a) self-reported osteoporosis instead of DXA category
fit <- try(svyglm(as.formula(paste("cvd ~ self_op +",
                                   paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + "))),
                  design = des, family = quasibinomial()), silent = TRUE)
if (!inherits(fit, "try-error")) {
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  k <- grep("^self_op", rownames(s))[1]
  if (!is.na(k)) {
    sens[[1]] <- data.frame(analysis = "self-reported osteoporosis", exposure = rownames(s)[k],
                            OR = exp(s[k, 1]), lcl = exp(ci[k, 1]), ucl = exp(ci[k, 2]),
                            p = s[k, 4], n = nrow(fit$model))
  }
}
# (b) sex-neutral T-score reference
fit <- try(svyglm(as.formula(paste("cvd ~ bmd_cat_all +",
                                   paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + "))),
                  design = des, family = quasibinomial()), silent = TRUE)
if (!inherits(fit, "try-error")) {
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  keep <- grep("^bmd_cat_all", rownames(s))
  sens[[length(sens) + 1]] <- data.frame(analysis = "sex-neutral T-score",
                                         exposure = rownames(s)[keep],
                                         OR = exp(s[keep, 1]), lcl = exp(ci[keep, 1]),
                                         ucl = exp(ci[keep, 2]), p = s[keep, 4],
                                         n = nrow(fit$model))
}
# (c) lumbar spine BMD in place of femoral neck (OR per 0.1 g/cm2 lower)
des_sp <- update(des, ls_low = -ls_bmd / 0.1)
fit <- try(svyglm(as.formula(paste("cvd ~ ls_low +",
                                   paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + "))),
                  design = des_sp, family = quasibinomial()), silent = TRUE)
if (!inherits(fit, "try-error")) {
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  sens[[length(sens) + 1]] <- data.frame(analysis = "lumbar spine BMD (per 0.1 g/cm2 lower)",
                                         exposure = "ls_low", OR = exp(s["ls_low", 1]),
                                         lcl = exp(ci["ls_low", 1]), ucl = exp(ci["ls_low", 2]),
                                         p = s["ls_low", 4], n = nrow(fit$model))
}
# (d) exclude participants with prevalent CVD for the mortality models is handled in
#     Table 4 M3 (prevalent CVD is a covariate); additionally report event counts
events <- data.frame(
  analysis = c("all-cause death", "CVD death", "prevalent CVD", "prevalent CHD", "prevalent stroke"),
  n_events = c(sum(d$dead, na.rm = TRUE), sum(d$cvd_death_num, na.rm = TRUE),
               sum(d$cvd, na.rm = TRUE), sum(d$chd, na.rm = TRUE), sum(d$stroke, na.rm = TRUE))
)
write.csv(do.call(rbind, sens), file.path(tabdir, "table5_sensitivity.csv"), row.names = FALSE)
write.csv(events, file.path(tabdir, "table6_event_counts.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Table 7: interaction tests
inter_p <- function(design, outcome, modvar) {
  form <- as.formula(paste(outcome, "~ bmd_cat *", modvar, "+",
                           paste(c("age", "sex", "race", "cycle", cov_metab), collapse = " + ")))
  fit <- try(svyglm(form, design = design, family = quasibinomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  tt <- try(regTermTest(fit, as.formula(paste0("~bmd_cat:", modvar))), silent = TRUE)
  if (inherits(tt, "try-error")) return(NA_real_)
  tt$p
}
inter <- data.frame(
  outcome = "cvd",
  interaction = c("bmd_cat x sex", "bmd_cat x age group"),
  p = c(inter_p(des, "cvd", "sex"), inter_p(des, "cvd", "age65"))
)
write.csv(inter, file.path(tabdir, "table7_interaction_tests.csv"), row.names = FALSE)
cat("\n=== interaction tests (prevalent CVD) ===\n"); print(inter, row.names = FALSE)

# ------------------------------------------- Table 8 and 9 live in dedicated scripts
# 04_rcs.R  -> table8_rcs_cvd_by_t_score.csv
# 05_km.R   -> table9_km_allcause_by_bmd.csv  (weighted Kaplan-Meier; svykm is
#              prohibitively slow here, so the curve is computed directly from the
#              survey weights and reported as a descriptive counterpart to Table 4)

# ---------------------------------------------------------------- sample flow
flow <- data.frame(
  step = c("aged >=40 in the 5 DXA cycles", "with femoral neck BMD",
           "and non-missing CVD status", "and mortality-eligible (final analytic n)"),
  n = c(18913, sum(d$fn_bmd > 0), sum(d$fn_bmd > 0 & !is.na(d$cvd)), nrow(d))
)
write.csv(flow, file.path(tabdir, "table0_sample_flow.csv"), row.names = FALSE)

cat("\n=== Table 3 (prevalent CVD, BMD category, M4) ===\n")
print(t3[t3$model == "M4 +mineral metabolism" & t3$subgroup == "All" & t3$exposure == "bmd_cat",
        c("term", "OR", "lcl", "ucl", "p", "n")], row.names = FALSE)
cat("\n=== Table 4 (mortality) ===\n")
print(t4, row.names = FALSE, digits = 3)
cat("\n=== Table 5 (sensitivity) ===\n")
print(do.call(rbind, sens), row.names = FALSE, digits = 3)
cat("\n=== weighted CVD prevalence by BMD category ===\n")
print(prev_out, row.names = FALSE, digits = 3)

