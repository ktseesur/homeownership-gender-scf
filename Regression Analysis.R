#Khongorzul Tseesuren
#Probit estimates for the paper "Who Holds the Keys? Gender and Homeownership Among Single Adults in the United States"
#Date: 11-23-2026
#Written for RStudio
#This file estimates the probit models and creates the regression table.
#It uses the analysis data created by "Data Construction.R" (2022 Survey of Consumer Finances).

#Remove all objects from the global environment
rm(list = ls())

# Set working directory (should be set to the main project folder)
# setwd("path/to/project-folder")  # the folder that contains Original Data/, Analysis Data/ and Documents/

# Install the following packages prior to using this script file:
if(!require("sandwich")) install.packages("sandwich")
if(!require("lmtest")) install.packages("lmtest")
if(!require("car")) install.packages("car")
if(!require("haven")) install.packages("haven")

# Load the following packages
library("sandwich")  # robust standard errors
library("lmtest")    # coefficient tests
library("car")       # variance inflation factors
library("haven")     # for reading the Stata file with the replicate weights

# Load analysis data file (created by Data Construction.R)
load("Analysis Data/scfdata.RData")


#### 1. Check the data ####

# Each family appears 5 times (5 implicates), so rows / 5 = families
nrow(scfdata) / 5
table(scfdata$implicate)

# Two samples: all unpartnered adults, and never-married adults only
all_data <- scfdata
nm_data  <- scfdata[scfdata$never_married == 1, ]
nrow(all_data) / 5   # families, all unpartnered
nrow(nm_data) / 5    # families, never married

# SCF replicate weights (file scf2022rw1s.zip from the Federal Reserve, unzipped
# into the Original Data folder). The file has one row per family (YY1) and
# 999 replicates. Each replicate is one more possible sample of families.
# WT1B1 to WT1B999 are the replicate weights and MM1 to MM999 count how many
# times the family was picked in that replicate.
rw_file <- list.files("Original Data", pattern = "rw1.*\\.dta$",
                      ignore.case = TRUE, full.names = TRUE)[1]
rw <- as.data.frame(read_dta(rw_file))
names(rw) <- tolower(names(rw))
rw[is.na(rw)] <- 0                      # a family that was not picked has weight 0

# Replicate weight = weight x number of times picked (one column per replicate)
rep_weights <- as.matrix(rw[, paste0("wt1b", 1:999)]) * as.matrix(rw[, paste0("mm", 1:999)])
rownames(rep_weights) <- rw$yy1

# Check: every family in the analysis data has replicate weights (should be TRUE)
all(scfdata$YY1 %in% rw$yy1)


#### 2. Models (controls added one step at a time) ####

# Column 1: only Female (the raw gap, as a probit)
f1 <- owner ~ female

# Column 2: add age, age squared, college, nonwhite, never married
f2 <- owner ~ female + age + age_sq + college + nonwhite + never_married

# Column 3: add children and Female x Children (tests H2)
f3 <- owner ~ female + age + age_sq + college + nonwhite + never_married +
  kids + female_kids

# Column 4: add log income, non-housing debt to income, credit denied
# This is the full model (preferred model)
f4 <- owner ~ female + age + age_sq + college + nonwhite + never_married +
  kids + female_kids + ln_income + nh_dti + credit_denied

# Column 5: the full model for never-married adults only
# (Never married is the same for everyone in this sample, so it is left out)
f5 <- owner ~ female + age + age_sq + college + nonwhite +
  kids + female_kids + ln_income + nh_dti + credit_denied


#### 3. Functions ####

# The SCF fills in missing answers five times (five implicates).
# So each model is estimated five times, once on each implicate, and combined.
# Standard errors follow the method in the SCF codebook:
#   coefficient       = average of the five coefficients
#   sampling variance = how much the coefficients move across the 999 replicate
#                       weights (estimated on the first implicate)
#   imputation variance = how much the coefficients move between the five implicates
#   total variance    = sampling variance + (6/5) x imputation variance

# Estimate one model with one set of weights
fit_one <- function(formula, d, type) {
  if (type == "probit") {
    glm(formula, family = quasibinomial(link = "probit"), data = d, weights = w,
        control = glm.control(maxit = 100))
  } else {
    lm(formula, data = d, weights = w)
  }
}

# Estimate one model on each implicate with the SCF weight (WGT).
# type = "probit" for the probit model, "ols" for the linear probability model.
# Robust (HC1) variances are also saved, only to compare in the robustness section.
fit_five <- function(formula, data, type = "probit") {
  fits <- list()
  for (i in 1:5) {
    d <- data[data$implicate == i, ]
    d$w <- d$WGT / mean(d$WGT)          # weight with mean 1, so weights add up to the number of families
    m <- fit_one(formula, d, type)
    fits[[i]] <- list(model = m, data = d, coef = coef(m),
                      robust = vcovHC(m, type = "HC1"))
  }
  fits
}

# Sampling variance from the 999 replicate weights (first implicate).
# The model is estimated 999 times, each time with a different replicate weight.
replicate_variance <- function(formula, fits, type = "probit") {
  d <- fits[[1]]$data
  rep_coefs <- sapply(1:999, function(r) {
    d$w <- rep_weights[as.character(d$YY1), r]
    dr  <- d[d$w > 0, ]                 # families picked in this replicate
    dr$w <- dr$w / mean(dr$w)
    coef(suppressWarnings(fit_one(formula, dr, type)))
  })
  cov(t(rep_coefs))                     # variance of the 999 sets of coefficients
}

# Combine the five implicates and the replicate-weight variance
combine_five <- function(formula, fits, type = "probit") {
  m        <- length(fits)
  coefs    <- sapply(fits, function(f) f$coef)            # one column per implicate
  q        <- rowMeans(coefs)                              # combined coefficient
  sampling <- replicate_variance(formula, fits, type)
  between  <- cov(t(coefs))                                # imputation variance
  total    <- sampling + (1 + 1 / m) * between             # combined variance
  se       <- sqrt(diag(total))
  z        <- q / se
  p        <- 2 * pnorm(-abs(z))
  # Robust standard errors with the same combining rule, for comparison only
  robust   <- Reduce("+", lapply(fits, function(f) f$robust)) / m + (1 + 1 / m) * between
  list(coef = q, se = se, z = z, p = p, vcov = total, robust_se = sqrt(diag(robust)))
}

# McFadden's pseudo R-squared (average of the five implicates)
pseudo_r2 <- function(fits) {
  mean(sapply(fits, function(f) 1 - f$model$deviance / f$model$null.deviance))
}

# Joint test: are all slope coefficients zero? (Wald chi-squared test)
joint_test <- function(comb) {
  b <- comb$coef[-1]
  V <- comb$vcov[-1, -1, drop = FALSE]
  chi2 <- as.numeric(t(b) %*% solve(V) %*% b)
  df <- length(b)
  c(chi2 = chi2, df = df, p = 1 - pchisq(chi2, df))
}

# Results table for one model (coefficient, standard error, z, p-value)
show_results <- function(comb) {
  round(cbind(Coefficient = comb$coef, Std.Error = comb$se,
              z = comb$z, p.value = comb$p), 4)
}


#### 4. Estimate the probit models ####

# Each line estimates the model 5 times (implicates) and 999 times (replicates),
# so this section takes a few minutes.
fits1 <- fit_five(f1, all_data); res1 <- combine_five(f1, fits1)
fits2 <- fit_five(f2, all_data); res2 <- combine_five(f2, fits2)
fits3 <- fit_five(f3, all_data); res3 <- combine_five(f3, fits3)
fits4 <- fit_five(f4, all_data); res4 <- combine_five(f4, fits4)
fits5 <- fit_five(f5, nm_data);  res5 <- combine_five(f5, fits5)

show_results(res1)
show_results(res2)
show_results(res3)
show_results(res4)   # preferred model, all unpartnered adults
show_results(res5)   # never-married adults only

# Model fit: McFadden's pseudo R-squared
fit_stats <- c(col1 = pseudo_r2(fits1), col2 = pseudo_r2(fits2), col3 = pseudo_r2(fits3),
               col4 = pseudo_r2(fits4), col5 = pseudo_r2(fits5))
round(fit_stats, 3)

# Joint significance of all explanatory variables (Wald test)
joint4 <- joint_test(res4)
joint5 <- joint_test(res5)
round(joint4, 4)
round(joint5, 4)


#### 5. Average partial effects and predicted probabilities ####

# Average partial effects = coefficient x scale factor.
# The scale factor is the (weighted) average of the normal density at each family's index.
ape_factor <- function(fits) {
  mean(sapply(fits, function(f) {
    xb <- predict(f$model)                 # linear index
    sum(f$data$w * dnorm(xb)) / sum(f$data$w)
  }))
}
ape4 <- res4$coef * ape_factor(fits4)
ape5 <- res5$coef * ape_factor(fits5)
round(cbind(APE_all_unpartnered = ape4), 4)
round(cbind(APE_never_married = ape5), 4)

# Predicted probability of owning for a woman and a man with the same controls.
# Every family keeps its own age, income, and other controls. Only Female and
# Children are changed. The gap is the woman's probability minus the man's.
predicted_gap <- function(fits, children) {
  one <- sapply(fits, function(f) {
    d_w <- f$data; d_w$female <- 1; d_w$kids <- children; d_w$female_kids <- children
    d_m <- f$data; d_m$female <- 0; d_m$kids <- children; d_m$female_kids <- 0
    p_w <- predict(f$model, newdata = d_w, type = "response")
    p_m <- predict(f$model, newdata = d_m, type = "response")
    c(woman = sum(f$data$w * p_w) / sum(f$data$w),
      man   = sum(f$data$w * p_m) / sum(f$data$w))
  })
  out <- rowMeans(one)
  c(out, gap = unname(out["woman"] - out["man"]))
}

# All unpartnered adults (preferred model): no children, then one child
gap4_nokids <- predicted_gap(fits4, children = 0)
gap4_onekid <- predicted_gap(fits4, children = 1)
round(rbind(no_children = gap4_nokids, one_child = gap4_onekid), 4)

# Never-married adults: no children, then one child
gap5_nokids <- predicted_gap(fits5, children = 0)
gap5_onekid <- predicted_gap(fits5, children = 1)
round(rbind(no_children = gap5_nokids, one_child = gap5_onekid), 4)


#### 6. Diagnostic checks and robustness ####

# (a) Multicollinearity: variance inflation factors (first implicate, full model)
# Age and age squared are related by construction, so their VIFs are high.
vif_model <- lm(f4, data = fits4[[1]]$data, weights = w)
vif_all <- vif(vif_model)
round(vif_all, 2)

# VIFs without age squared and the interaction (to see the other variables)
vif_model2 <- lm(owner ~ female + age + college + nonwhite + never_married +
                   kids + ln_income + nh_dti + credit_denied,
                 data = fits4[[1]]$data, weights = w)
vif_simple <- vif(vif_model2)
round(vif_simple, 2)

# (b) Robustness: linear probability model (OLS), same weights and standard errors
lpm4_fits <- fit_five(f4, all_data, type = "ols"); lpm4 <- combine_five(f4, lpm4_fits, type = "ols")
lpm5_fits <- fit_five(f5, nm_data,  type = "ols"); lpm5 <- combine_five(f5, lpm5_fits, type = "ols")
show_results(lpm4)
show_results(lpm5)

# Compare the Female coefficient: probit average partial effect and OLS
round(c(probit_APE_all = unname(ape4["female"]), OLS_all = unname(lpm4$coef["female"]),
        probit_APE_never_married = unname(ape5["female"]),
        OLS_never_married = unname(lpm5$coef["female"])), 4)

# (c) Robustness: full model without families with income of $1 million or more
# A family is dropped if its income is $1 million or more in any of the five implicates
rich_families <- unique(all_data$YY1[all_data$ln_income >= log(1e6)])
length(rich_families)                 # families dropped
no_rich <- all_data[!(all_data$YY1 %in% rich_families), ]
nrow(no_rich) / 5                     # families kept
fits4b <- fit_five(f4, no_rich); res4b <- combine_five(f4, fits4b)
show_results(res4b)

# (d) Robustness: replicate-weight standard errors next to robust (HC1) standard errors
# for the preferred model. The coefficients are the same; only the standard errors differ.
se_compare <- cbind(Coefficient = res4$coef, SE_replicate = res4$se, SE_robust = res4$robust_se,
                    p_replicate = res4$p,
                    p_robust = 2 * pnorm(-abs(res4$coef / res4$robust_se)))
round(se_compare, 4)


#### 7. Regression table ####

# The table is written as an HTML file with base R (no extra package needed).
# Each cell shows the combined coefficient with stars, and the standard error
# (from the replicate weights) in parentheses below it.

results_list <- list(res1, res2, res3, res4, res5, lpm4)
column_names <- c("(1) All", "(2) All", "(3) All", "(4) All",
                  "(5) Never married", "(6) All, OLS")
families     <- c(rep(nrow(all_data) / 5, 4), nrow(nm_data) / 5, nrow(all_data) / 5)
pseudo_row   <- c(sprintf("%.3f", fit_stats), "")

# Variable names in the models and the meaningful names shown in the table
var_names  <- c("female", "age", "age_sq", "college", "nonwhite", "never_married",
                "kids", "female_kids", "ln_income", "nh_dti", "credit_denied",
                "(Intercept)")
var_labels <- c("Female", "Age", "Age squared", "College", "Nonwhite", "Never married",
                "Children", "Female x Children", "Log income",
                "Non-housing debt to income", "Credit denied", "Constant")

# Significance stars from the p-value
stars <- function(p) {
  if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else ""
}

# Three decimals, or more when the number is very small (age squared)
fmt <- function(x) {
  if (abs(x) < 0.0001) sprintf("%.5f", x) else if (abs(x) < 0.001) sprintf("%.4f", x) else sprintf("%.3f", x)
}

# One table cell: coefficient with stars, standard error below (empty if not in the model)
make_cell <- function(res, v) {
  if (!(v %in% names(res$coef))) return("")
  paste0(fmt(res$coef[v]), stars(res$p[v]), "<br>(", fmt(res$se[v]), ")")
}

# One table row from a label and its cells
make_row <- function(label, cells) {
  paste0("<tr><td style='text-align:left'>", label, "</td>",
         paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
}

table_rows <- character(0)
for (k in seq_along(var_names)) {
  cells <- sapply(results_list, make_cell, v = var_names[k])
  table_rows <- c(table_rows, make_row(var_labels[k], cells))
}

html <- c(
  "<html><head><style>td {padding: 1px 12px}</style></head>",
  "<body style='font-family: Times New Roman; font-size: 12pt'>",
  "<table style='border-collapse: collapse; text-align: center'>",
  "<caption style='text-align:left'><b>Table 4: Probit estimates of homeownership</b><br>",
  "Dependent variable: Homeowner (1 = owns primary residence)</caption>",
  "<tr style='border-top: 2px solid black; border-bottom: 1px solid black'><td></td>",
  paste0("<td>", column_names, "</td>", collapse = ""), "</tr>",
  table_rows,
  paste0("<tr style='border-top: 1px solid black'><td style='text-align:left'>Families</td>",
         paste0("<td>", format(families, big.mark = ","), "</td>", collapse = ""), "</tr>"),
  paste0("<tr style='border-bottom: 2px solid black'><td style='text-align:left'>Pseudo R-squared</td>",
         paste0("<td>", pseudo_row, "</td>", collapse = ""), "</tr>"),
  "</table>",
  "<p>Notes: Columns 1 to 5 are probit coefficients. Column 6 is a linear probability model (OLS).<br>",
  "Standard errors, in parentheses, use the 999 SCF replicate weights and the five implicates.<br>",
  "Estimates use the SCF weight (WGT) and combine the five implicates.<br>",
  "* p &lt; 0.10, ** p &lt; 0.05, *** p &lt; 0.01</p>",
  "</body></html>")

writeLines(html, "Documents/regression_table.html")


#### 8. Save all results in one text file ####

sink("Documents/Regression_output.txt")
cat("Probit results, 2022 SCF, unpartnered adults\n")
cat("Families, all unpartnered:", nrow(all_data) / 5, "\n")
cat("Families, never married:", nrow(nm_data) / 5, "\n\n")
cat("Column 1 (Female only)\n");                 print(show_results(res1))
cat("\nColumn 2 (+ age, education, race, never married)\n"); print(show_results(res2))
cat("\nColumn 3 (+ children, Female x Children)\n");          print(show_results(res3))
cat("\nColumn 4 (+ income, debt, credit denied; preferred)\n"); print(show_results(res4))
cat("\nColumn 5 (never-married adults, full model)\n");        print(show_results(res5))
cat("\nPseudo R-squared\n");                     print(round(fit_stats, 3))
cat("\nJoint Wald test, column 4\n");            print(round(joint4, 4))
cat("\nJoint Wald test, column 5\n");            print(round(joint5, 4))
cat("\nAverage partial effects, column 4\n");    print(round(ape4, 4))
cat("\nAverage partial effects, column 5\n");    print(round(ape5, 4))
cat("\nPredicted probabilities, all unpartnered (woman, man, gap)\n")
print(round(rbind(no_children = gap4_nokids, one_child = gap4_onekid), 4))
cat("\nPredicted probabilities, never married (woman, man, gap)\n")
print(round(rbind(no_children = gap5_nokids, one_child = gap5_onekid), 4))
cat("\nVariance inflation factors, full model\n");            print(round(vif_all, 2))
cat("\nVariance inflation factors, without age squared and interaction\n"); print(round(vif_simple, 2))
cat("\nLinear probability model (OLS), all unpartnered\n");   print(show_results(lpm4))
cat("\nLinear probability model (OLS), never married\n");     print(show_results(lpm5))
cat("\nFull probit without incomes of $1 million or more (families:", nrow(no_rich) / 5, ")\n")
print(show_results(res4b))
cat("\nPreferred model: replicate-weight and robust standard errors\n")
print(round(se_compare, 4))
sink()
