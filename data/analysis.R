load("data/MIexample.RData")

# fileConn <- file("E:/GitHub/EpiShinyLiveM/data/analyticData2.txt")
# dput(dat.full, fileConn)
# close(fileConn)



# Load required packages
library(dplyr)
library(kableExtra)
library(tableone)
library(survey)
library(Publish)
library(DataExplorer)
library(mice)
library(mitools)
dim(dat.full)
head(dat.full)
# Drop < 20 years
dat.with.miss <- subset(dat.full, age >= 20)

# Frequency for outcome and exposure 
table(dat.with.miss$cvd, useNA = "always") # 6 missing
table(dat.with.miss$rheumatoid, useNA = "always") # 1375 missing

# Drop missing in outcome and exposure - dataset with missing values only in covariates
dat.analytic <- dat.with.miss[complete.cases(dat.with.miss$cvd),]
dat.analytic <- dat.analytic[complete.cases(dat.analytic$rheumatoid),]
nrow(dat.analytic)
dat.full$ineligible <- 1
dat.full$ineligible[dat.full$studyid %in% dat.analytic$studyid] <- 0
table(dat.full$ineligible, useNA = "always")
# Categorical age
dat.analytic$age.cat <- with(dat.analytic, ifelse(age >= 20 & age < 50, "20-49", 
                                                  ifelse(age >= 50 & age < 65, "50-64", "65+")))
dat.analytic$age.cat <- factor(dat.analytic$age.cat, levels = c("20-49", "50-64", "65+"))
table(dat.analytic$age.cat, useNA = "always")

# Recode rheumatoid to arthritis
dat.analytic$arthritis <- car::recode(dat.analytic$rheumatoid, " 'No' = 'No arthritis';
                                      'Yes' = 'Rheumatoid arthritis' ", as.factor = T)
table(dat.analytic$arthritis, useNA = "always")

# Keep only relevant variables
vars <-  c("studyid", "survey.weight", "psu", "strata", "cvd", "arthritis", "age.cat", 
           "sex", "education", "race", "income", "bmi", "smoking", "htn", "diabetes")
dat.analytic2 <- dat.analytic[, vars]
# Create Table 1
vars <- c("arthritis", "age.cat", "sex", "education", "race", "income", "bmi", "smoking",
          "htn", "diabetes")
tab1 <- CreateTableOne(vars = vars, strata = "cvd", data = dat.analytic2, includeNA = F,
                       addOverall = T, test = F)
print(tab1, format = "f", showAllLevels = T)
DataExplorer::plot_missing(dat.analytic2)
# Step 0: Set imputation model
ini <- mice(data = dat.analytic2, maxit = 0, print = FALSE)
pred <- ini$pred

# Use the strata variable as an auxiliary variable in the imputation model
pred["strata",] <- 0

# Do not use survey weight or PSU variable as auxiliary variables
pred[,"studyid"] <- pred["studyid",] <- 0
pred[,"psu"] <- pred["psu",] <- 0
pred[,"survey.weight"] <- pred["survey.weight",] <- 0

# Set imputation method
meth <- ini$meth
meth["bmi"] <- "pmm"
meth["income"] <- "pmm"
meth
# Step 1: impute the incomplete data
imputation <- mice(data = dat.analytic2,
                   seed = 123,
                   predictorMatrix = pred,
                   method = meth,
                   m = 5,
                   maxit = 3,
                   print = FALSE)
impdata <- mice::complete(imputation, action="long")

table(impdata$age.cat)
#Remove .id variable from the model as it was created in an intermediate step
impdata$.id <- NULL

# Create an indicator of eligible subjects 
impdata$ineligible <- 0

# Number of subjects
nrow(impdata)
DataExplorer::plot_missing(impdata)
# Number of ineligible subjects
#dat.full$ineligible <- 1
#dat.full$ineligible[dat.full$studyid %in% dat.analytic$studyid] <- 0
table(dat.full$ineligible, useNA = "always")
# Subset for ineligible
dat.ineligible <- subset(dat.full, ineligible == 1)

# Create m = 5 datasets with .imp 1 to m = 5
dat31 <- dat.ineligible; dat31$.imp <- 1
dat32 <- dat.ineligible; dat32$.imp <- 2
dat33 <- dat.ineligible; dat33$.imp <- 3
dat34 <- dat.ineligible; dat34$.imp <- 4
dat35 <- dat.ineligible; dat35$.imp <- 5
# Stacked data for ineligible subjects
dat.ineligible.stacked <- rbind(dat31, dat32, dat33, dat34, dat35)
DataExplorer::plot_missing(dat.ineligible.stacked)
names(impdata)
names(dat.ineligible.stacked)
# Categorical age
summary(dat.ineligible.stacked$age)
dat.ineligible.stacked$age.cat <- with(dat.ineligible.stacked, 
                                       ifelse(age >= 20 & age < 50, "20-49", 
                                              ifelse(age >= 50 & age < 65, "50-64", 
                                                     ifelse(age >= 65, "65+", NA))))
dat.ineligible.stacked$age.cat <- factor(dat.ineligible.stacked$age.cat, 
                                         levels = c("20-49", "50-64", "65+"))
table(dat.ineligible.stacked$age.cat, useNA = "always")
# Recode arthritis
dat.ineligible.stacked$arthritis <- car::recode(dat.ineligible.stacked$rheumatoid, 
                                                " 'No' = 'No arthritis'; 'Yes' = 
                                                'Rheumatoid arthritis' ", as.factor = T)
# Variable names in the imputed dataset
vars <- names(impdata) 

# Set up the dataset for ineligible - same variables as impdata
dat.ineligible.stacked <- dat.ineligible.stacked[, vars]
impdata2 <- rbind(impdata, dat.ineligible.stacked)
impdata2 <- impdata2[order(impdata2$.imp, impdata2$studyid),]
dim(impdata2)
m <- 5
allImputations <- imputationList(lapply(1:m, function(n) subset(impdata2, subset=.imp==n)))

# Step 2: Survey data analysis
w.design0 <- svydesign(ids = ~psu, 
                       weights = ~survey.weight, 
                       strata = ~strata,
                       data = allImputations, 
                       nest = TRUE) # Design on full data
w.design <- subset(w.design0, ineligible == 0) # Subset the design
dim(w.design)
# Design-adjusted logistic regression
fit <- with(w.design, svyglm(I(cvd == "Yes") ~ arthritis + age.cat + sex + education + 
                               race + income + bmi + smoking + htn + diabetes, 
                             family = quasibinomial))
res <- exp(as.data.frame(cbind(coef(fit[[1]]),
                               coef(fit[[2]]),
                               coef(fit[[3]]),
                               coef(fit[[4]]),
                               coef(fit[[5]]))))
names(res) <- paste("OR from m =", 1:5)
round(t(res),2)
# Step 3: Pooled estimates
pooled.estimates <- MIcombine(fit)
Estimate <- (exp(pooled.estimates$coefficients))
OR <- as.data.frame(Estimate)
CI <- (exp(confint(pooled.estimates)))
OR <- cbind(OR, CI)
names(OR) <- c("Estimate", "lower", "upper")
OR <- rownames_to_column(OR, var = "term")
OR$term <- factor(OR$term, levels = OR$term)
OR

ggplot(OR, aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")


list.res <- vector(mode='list', length=5)

for (i in 1:5){
  df.resid.val <- degf(w.design[[1]][[i]])
  
  # Extract coefficients and confidence intervals
  coef_summary <- coef(summary(fit[[i]], df.resid = df.resid.val))
  ci_values <- confint.default(fit[[i]], df.resid = df.resid.val)
  
  # Combine and exponentiate
  combined <- cbind(coef_summary[, "Estimate"], ci_values)
  colnames(combined) <- c("Estimate", "lower", "upper")
  combined <- exp(combined)
  
  # Convert to data frame and add row names as column
  combined_df <- as.data.frame(combined)
  combined_df <- rownames_to_column(combined_df, var = "term")
  combined_df$term <- factor(combined_df$term, levels = combined_df$term)
  
  list.res[[i]] <- combined_df
}

list.res$MI <- OR

ggplot(list.res[[1]], aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")

ggplot(list.res[[2]], aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")

ggplot(list.res[[3]], aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")

ggplot(list.res[[4]], aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")

ggplot(list.res[[5]], aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")
ggplot(list.res$MI, aes(y = term, x = Estimate, xmin = lower, xmax = upper)) +
  geom_point() +  # for the point estimates
  geom_errorbarh(height = 0.2) +  # for the horizontal confidence intervals
  theme_minimal() +
  labs(x = "Odds Ratio (OR)", y = "", title = "Forest Plot of Regression Coefficients")


# edit(list.res)
