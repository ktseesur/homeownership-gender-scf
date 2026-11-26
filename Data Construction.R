#Khongorzul Tseesuren
#Data construction for the paper "Who Holds the Keys? Gender and Homeownership Among Single Adults in the United States"
#Date: 11-13-2026
#Analysis data from the 2022 Survey of Consumer Finances

#Remove all objects from the global environment
rm(list = ls())

# Set working directory (should be set to the main project folder)
# setwd("path/to/project-folder")  # the folder that contains Original Data/, Analysis Data/ and Documents/

# Install the following packages prior to using this script file:
if(!require("dplyr")) install.packages("dplyr")
if(!require("labelled")) install.packages("labelled")
if(!require("haven")) install.packages("haven")

# Load the following packages
library(dplyr) # filtering, merging and creating variables
library(labelled) # for variable labels
library(haven) # for reading the Stata file


#Load data

# 2022 SCF, Summary Extract Public Data (SCFP2022.csv)
scfdata <- read.csv("Original Data/SCFP2022.csv")
names(scfdata) <- toupper(names(scfdata)) # upper case to match the SCF codebook

# 2022 SCF, Full Public Data Set (p22i6.dta, unzipped from scf2022s.zip)
# Only the row ID and detailed marital status (X8023) are needed:
# 1 = married, 2 = partner, 3 = separated, 4 = divorced, 5 = widowed, 6 = never married
fulldata <- read_dta("Original Data/p22i6.dta", col_select = c(y1, x8023))
names(fulldata) <- toupper(names(fulldata))
fulldata$X8023 <- as.numeric(fulldata$X8023)


#Check the data

# Size of each file (each family appears 5 times, once per implicate)
dim(scfdata)
dim(fulldata)

# Variables used in this project and their types
vars <- c("YY1", "Y1", "WGT", "HOUSECL", "HHSEX", "MARRIED", "AGE", "INCOME",
          "EDCL", "KIDS", "RACECL4", "DEBT", "MRTHEL", "RESDBT", "TURNFEAR")
str(scfdata[, vars])

# Number of families and number of rows per family (should be 5)
length(unique(scfdata$YY1))
table(table(scfdata$YY1))

# Missing values (the SCF imputes missing answers, so there should be none)
colSums(is.na(scfdata[, vars]))
sum(is.na(fulldata$X8023))

# Duplicate observations (Y1 is the row ID)
sum(duplicated(scfdata$Y1))
sum(duplicated(fulldata$Y1))

# Coding of categorical variables
table(scfdata$HOUSECL)   # 1 = owns home, 2 = otherwise
table(scfdata$HHSEX)     # 1 = male, 2 = female
table(scfdata$MARRIED)   # 1 = married or living with partner, 2 = neither
table(scfdata$EDCL)      # 1 to 4 education groups
table(scfdata$RACECL4)   # 1 to 4 race and ethnicity groups
table(scfdata$TURNFEAR)  # 1 = turned down for credit or feared denial
table(fulldata$X8023)    # 1 to 6 marital status

# Ranges of numerical variables
summary(scfdata$AGE)
summary(scfdata$INCOME)
summary(scfdata$KIDS)
summary(scfdata$DEBT)
summary(scfdata$WGT)

# Weights sum to about the number of U.S. families (WGT is already divided by 5)
sum(scfdata$WGT)

# Weighted homeownership rate for all families (Fed publishes 66.1% for 2022)
100 * sum(scfdata$WGT * (scfdata$HOUSECL == 1)) / sum(scfdata$WGT)


#Merge marital status into the summary extract

# Both files share the row ID Y1
scfdata <- scfdata %>%
  left_join(fulldata, by = "Y1")

# Check the merge (same number of rows, no missing marital status)
nrow(scfdata)
sum(is.na(scfdata$X8023))

# MARRIED = 1 should only match X8023 = 1 or 2, and MARRIED = 2 only 3 to 6
table(MARRIED = scfdata$MARRIED, X8023 = scfdata$X8023) / 5


#Keep families without a partner (MARRIED = 2)
# The SCF lists the man as the reference person in mixed-sex couples,
# so the sex of the reference person only shows gender when there is no partner

scfdata <- scfdata %>%
  filter(MARRIED == 2)

length(unique(scfdata$YY1))   # families in the sample
table(scfdata$HHSEX) / 5      # families by sex


#Drop families with zero income (log income is undefined at zero)

sum(scfdata$INCOME == 0) / 5  # families dropped

scfdata <- scfdata %>%
  filter(INCOME > 0)

length(unique(scfdata$YY1))   # families in the final sample
nrow(scfdata)                 # rows (families x 5 implicates)


# Generate analysis variables
# Dollar amounts in the SCF are already in 2022 dollars, so no inflation adjustment is needed

scfdata <- scfdata %>%
  mutate(
    implicate     = Y1 - 10 * YY1,                     # implicate number 1 to 5
    owner         = as.integer(HOUSECL == 1),          # 1 = owns home
    female        = as.integer(HHSEX == 2),            # 1 = female reference person
    age           = AGE,
    age_sq        = AGE^2,
    age_group     = cut(AGE, breaks = c(-Inf, 34, 44, 54, 64, 74, Inf),
                        labels = c("Under 35", "35-44", "45-54",
                                   "55-64", "65-74", "75 and older")),
    income_k      = INCOME / 1000,                     # thousands of 2022 dollars (tables only)
    ln_income     = log(INCOME),                       # income measure in the model
    educ_group    = factor(EDCL, levels = 1:4,
                           labels = c("No high school diploma",
                                      "High school diploma or GED",
                                      "Some college or associate degree",
                                      "Bachelor's degree or higher")),
    college       = as.integer(EDCL == 4),
    kids          = KIDS,
    has_kids      = as.integer(KIDS > 0),
    female_kids   = female * kids,                     # female x number of children
    race_group    = factor(RACECL4, levels = 1:4,
                           labels = c("White non-Hispanic",
                                      "Black non-Hispanic",
                                      "Hispanic or Latino",
                                      "Other or multiple races")),
    nonwhite      = as.integer(RACECL4 != 1),
    nh_dti        = (DEBT - MRTHEL - RESDBT) / INCOME, # non-housing debt to income (debt not on real estate)
    credit_denied = as.integer(TURNFEAR == 1),         # 1 = turned down or feared denial
    marital       = factor(X8023, levels = 3:6,
                           labels = c("Separated", "Divorced",
                                      "Widowed", "Never married")),
    never_married = as.integer(X8023 == 6)
  )

# Check the new variables
summary(scfdata[, c("owner", "female", "age", "age_sq", "ln_income", "college",
                    "kids", "female_kids", "nonwhite", "nh_dti",
                    "credit_denied", "never_married")])

# Non-housing debt should not be negative
sum(scfdata$nh_dti < 0)

# Never-married families by sex, and how many have children
table(female = scfdata$female, never_married = scfdata$never_married) / 5
table(female = scfdata$female[scfdata$never_married == 1],
      has_kids = scfdata$has_kids[scfdata$never_married == 1]) / 5

# Very high incomes are real families (the SCF oversamples wealthy families),
# so they are kept
sum(scfdata$INCOME >= 1e6) / 5


#Label each variable

var_label(scfdata) <- list(
  YY1           = "SCF family ID",
  Y1            = "SCF row ID (family ID x 10 + implicate)",
  implicate     = "SCF implicate number (1 to 5)",
  WGT           = "SCF sample weight (already divided by 5)",
  owner         = "Homeowner (1 = owns primary residence)",
  female        = "Female reference person (1 = female)",
  age           = "Age of reference person (years)",
  age_sq        = "Age squared",
  age_group     = "Age group of reference person",
  income_k      = "Family income in 2021 (thousands of 2022 dollars)",
  ln_income     = "Log of family income",
  educ_group    = "Education of reference person",
  college       = "Bachelor's degree or higher (1 = yes)",
  kids          = "Number of children in the household",
  has_kids      = "Has children in the household (1 = yes)",
  female_kids   = "Female x number of children",
  race_group    = "Race and ethnicity of reference person",
  nonwhite      = "Nonwhite or Hispanic reference person (1 = yes)",
  nh_dti        = "Non-housing debt divided by family income",
  credit_denied = "Turned down for credit or did not apply for fear of denial (1 = yes)",
  marital       = "Marital status of reference person",
  never_married = "Never married (1 = yes)"
)


#Keep selected variables
# WGT is kept because all statistics and models use the SCF weights

scfdata <- scfdata %>%
  select(YY1, Y1, implicate, WGT, owner, female, age, age_sq, age_group,
         income_k, ln_income, educ_group, college, kids, has_kids, female_kids,
         race_group, nonwhite, nh_dti, credit_denied, marital,
         never_married)

# No missing values in the final data
colSums(is.na(scfdata))

#Save R data file in Analysis Data folder
save(scfdata, file = "Analysis Data/scfdata.RData")


#Notes on the final data
# Only 18 never-married men have children in the household (143 never-married
#   women do), so all unpartnered families are kept to test children in the
#   larger sample too.
# Outliers in debt and income are kept because they are real families; WGT
#   accounts for the oversampling of wealthy families.
# The SCF has no credit scores, so credit_denied is used instead.
# Wealth is not included because home equity is part of wealth.
# The SCF has no measure of job mobility.
