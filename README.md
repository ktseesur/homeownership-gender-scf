# Who Holds the Keys? Gender and Homeownership Among Single Adults in the United States

R code for my paper (preliminary draft, November 2026), which asks whether unpartnered women are less likely than unpartnered men to own a home, and whether the gap differs when children live in the household. The data are the 2022 Survey of Consumer Finances (SCF).

Author: Khongorzul Tseesuren, Eastern Michigan University

## Scripts

1. `Data Construction.R` reads the SCF Summary Extract and the detailed marital status (X8023) from the Full Public Data Set, merges them by row ID, checks the data (row counts, missing values, duplicates, weight totals, the weighted homeownership rate), keeps families without a partner and with positive income, creates the analysis variables, and saves `Analysis Data/scfdata.RData`.
2. `Regression Analysis.R` estimates the probit models (Table 4 of the paper) on each of the five implicates with the SCF weight, computes standard errors from the 999 SCF replicate weights combined with the imputation variance, and computes average partial effects and predicted probabilities (Table 5, Panel A). It also runs the robustness checks (linear probability model, without incomes of $1 million or more, robust standard errors) and writes the regression table and a text file of all results.

## Data (not included in this repository)

Download these from the Federal Reserve's 2022 SCF page (https://www.federalreserve.gov/econres/scf_2022.htm) and put them in `Original Data/`:

- `SCFP2022.csv`: Summary Extract Public Data
- `p22i6.dta`: Full Public Data Set (unzipped from `scf2022s.zip`)
- the replicate weight file (unzipped from `scf2022rw1s.zip`)

## Folder layout

```
project-folder/
  Data Construction.R
  Regression Analysis.R
  Original Data/     the three SCF files
  Analysis Data/     created data (scfdata.RData)
  Documents/         regression_table.html and Regression_output.txt
```

Set the working directory to `project-folder/` before running.

## How to run

Run `Data Construction.R` first, then `Regression Analysis.R`. The second script takes a few minutes, because each model is estimated 999 times with the replicate weights.

## R packages

dplyr, labelled, haven, sandwich, lmtest, car
