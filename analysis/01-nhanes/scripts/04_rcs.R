# Step 4 - restricted cubic spline: OR for prevalent CVD across the sex-specific
# femoral neck T-score, referenced to T = 0. Contrasts use the fitted model's
# design-based covariance, so the confidence band inherits the survey design.
#
# Input : data/processed/nhanes/analysis_dataset.csv
# Output: results/tables/01-nhanes/table8_rcs_cvd_by_t_score.csv
suppressPackageStartupMessages({
  library(survey)
  library(splines)
})

tabdir <- "results/tables/01-nhanes"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv("data/processed/nhanes/analysis_dataset.csv", stringsAsFactors = FALSE)

d$bmd_cat <- factor(d$bmd_cat, levels = c("Normal", "Osteopenia", "Osteoporosis"))
d$sex <- factor(d$sex, levels = c("Male", "Female"))
d$race <- factor(d$race); d$race <- relevel(d$race, ref = "NH White")
d$smoking <- factor(d$smoking, levels = c("Never", "Former", "Current"))
d$education <- factor(d$education, levels = c(">High school", "High school/GED", "<High school"))
d$cycle <- factor(d$cycle)
d$strata2 <- factor(paste(d$cycle, d$strata, sep = "_"))
d$wt5 <- d$wt_mec / 5
d$cvd <- as.numeric(d$cvd)
d$egfr <- as.numeric(d$egfr)

des <- svydesign(ids = ~psu, strata = ~strata2, weights = ~wt5, data = d, nest = TRUE)

cov_all <- c("age", "race", "education", "pir", "bmi", "smoking",
             "diabetes", "hypertension", "total_chol", "hdl", "egfr",
             "albumin", "calcium", "phosphate", "alp")
form <- as.formula(paste0("cvd ~ ns(t_score, df = 3) + ",
                          paste(c("age", "sex", "race", "cycle", cov_all), collapse = " + ")))
fit <- svyglm(form, design = des, family = quasibinomial())
cat("model fitted, n =", nrow(fit$model), "\n")

grid <- seq(-4.5, 2.5, by = 0.1)
nd <- des$variables[rep(1, length(grid)), , drop = FALSE]
for (v in c("age", "pir", "bmi", "total_chol", "hdl", "egfr", "albumin", "calcium",
            "phosphate", "alp")) nd[[v]] <- median(des$variables[[v]], na.rm = TRUE)
for (v in c("sex", "race", "education", "smoking", "cycle")) {
  best <- names(which.max(table(des$variables[[v]])))
  nd[[v]] <- factor(best, levels = levels(des$variables[[v]]))
}
nd$hypertension <- 0; nd$diabetes <- 0
nd$t_score <- grid

ref <- nd[1, , drop = FALSE]; ref$t_score <- 0
X <- model.matrix(delete.response(terms(fit)), nd)
Xr <- model.matrix(delete.response(terms(fit)), ref)
beta <- coef(fit)
common <- intersect(colnames(X), names(beta))
X <- X[, common, drop = FALSE]; Xr <- Xr[, common, drop = FALSE]
miss <- setdiff(names(beta), common)
if (length(miss) > 0) {
  addcol <- matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss))
  X <- cbind(X, addcol)
}
X <- X[, names(beta), drop = FALSE]; Xr <- Xr[, names(beta), drop = FALSE]
D <- sweep(X, 2, Xr[1, ])
lp <- as.vector(D %*% beta)
V <- vcov(fit)
se <- sqrt(rowSums((D %*% V) * D))
out <- data.frame(t_score = grid, OR = exp(lp),
                  lcl = exp(lp - 1.96 * se), ucl = exp(lp + 1.96 * se))
# reference-anchored p for non-linearity: Wald test on spline terms
wt <- regTermTest(fit, ~ns(t_score, df = 3))
cat("spline term p =", format.pval(wt$p, digits = 3), "\n")
write.csv(out, file.path(tabdir, "table8_rcs_cvd_by_t_score.csv"), row.names = FALSE)
cat("->", file.path(tabdir, "table8_rcs_cvd_by_t_score.csv"), "\n")
print(head(out, 3)); print(tail(out, 3))

