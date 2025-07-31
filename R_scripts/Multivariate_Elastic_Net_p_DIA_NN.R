
## Script generated on 16/03/23
## Last update on 05/06/25
## Purpose: Identification of a subset of proteins (biomarker candidates) which
## are associated with both the disease causing CTG repeat expansion and a 
## gold standard clinical phenotype assessment (6MWT)

###############
## Libraries ##
###############

library(glmnet)  #V 4.1.4
library(mice)    #V 3.14.0
library(boot)    #V 1.3.28.1
library(dplyr)   #V 1.0.8
library(Metrics) #V 0.1.4
library(readxl)  #V 1.4.0
library(MASS)    #V 7.3.53
library(gvlma)   #V 1.0.0.3


###############
## Functions ##
###############

## Function to run cv.glmnet on bootstrapped datasets and return coefficients
fitglmnet <- function(data, indices, alpha, nfolds, dvariables, family){
  
  # Allow boot to select a bootstrapped dataset
  d <- data[indices,]
  
  # Matrix conversion & splitting dependent variable
  x <- as.matrix(d[,-which(names(d) %in% dvariables)])
  y <- as.matrix(d[,dvariables], drop=F)

  # Fit cross validated model
  cvfit <- cv.glmnet(x, y, type.measure = "mse", family=family,      
                     alpha=alpha, nfolds=nfolds, 
                     standardize = TRUE, standardize.response=TRUE)
  
  # Obtain & return coefficients
  coef_m <- coef(cvfit, s="lambda.1se")
  
  # In a multivariate setting, the same predictors are selected for all dependent 
  # variables. Hence, for the VIP calculating it suffices to obtain one set.
  if (length(dvariables) > 1){
    coef_m <- coef_m[[1]]                                                       
  }
  # Remove intercept
  coef_m <- coef_m[-1,]                                                         
  return(coef_m)
}

## Function to apply boxcox transformation with optimal lambda
boxcox_transformation <- function(x){
  lm_b <- lm(x ~ 1, y=TRUE)
  model <- boxcox(lm_b, plotit = FALSE)
  lambda <- model$x[which.max(model$y)]
  if (lambda == 0){
    new_x <- log(x)
  }
  else{
    new_x <- (((x^lambda)-1) / lambda)
  }
  return(new_x)
}

## Function to calculate (out of sample) R-squared values
rsq <- function(observed, predicted){
  # p = number of predictors excluding constant
  # N = total sample size
  
  # R-sq
  SSRes <- sum((observed - predicted)^2)
  SSTot <- sum((observed - mean(observed))^2)
  Rsq <- round(1 - (SSRes / SSTot), 4)

  return(Rsq)
}


###################
## Load datasets ##
###################

## External proteomics validation dataset
load(file = "all_p_ext_val_res.RDATA")
cdf_e <- cdf
ext_p_ids <- colnames(cdf_e)[8:ncol(cdf_e)]
rm(c, cdf, CTG_res_df, SMWT_res_df)

# Filtered, normalized and imputed DIA-NN protein quantification with namekey
load(file = "normalized_protein_intensities_p_DIA_NN.RDATA")
rm(rpc)

# Results of statistical analyses
load(file = "DIA-NN_p_Protein_statistics.RDATA")


##################################################################
## Multivariate modeling: CTG repeat & SMWT predictor selection ##
##################################################################

set.seed(1)

## Select significant predictors & assign dependent variables
SMWT <- outcome_coef[["SMWT"]]
SMWT <- SMWT[SMWT$SMWT_FDR < 0.05,]

CTG <- Mode_res_df[Mode_res_df$Mode_FDR < 0.05,]
table(rownames(SMWT) %in% rownames(CTG))
predictors <- intersect(rownames(SMWT), rownames(CTG))
dvariables <- c("SMWT", "Mode")

## Filter predictors for unambiguous identifications in external dataset
predictors <- predictors[predictors %in% ext_p_ids]

# Change protein names according to namekey [generated in protein_analysis.R]
prot_pred <- c()
for (prot in predictors){
  prot_pred <- c(prot_pred, namekey$Permuted[namekey$Original==prot])
}

## Generate dataframe to store selected predictors per bootstrapped distribution
selected_pred <- data.frame(matrix(ncol=length(prot_pred), nrow=0))
colnames(selected_pred) <- prot_pred

## For each imputed dataset, fit bootstrap enhanced Elastic-Net models
for (i in 1:10){
  print(paste("Currently fitting imputed dataset", i))
  
  # Obtain relevant dataset and filter for relevant predictors
  idf <- complete(rpc_Ti, i)
  idf <- idf[, colnames(idf) %in% prot_pred] 
  
  # Merge protein data with relevant phenotype data using patient ID's
  table(rownames(idf) == rownames(pdf))
  idf <- merge(idf, pdf[,c("SMWT", "Mode", "Visit")], by="row.names")
  rownames(idf) <- idf$Row.names
  idf$Row.names <- NULL
  
  # Remove cases with missing dependent variable information                    ## Removes 15 cases
  idf <- idf[!(is.na(idf$Mode) | is.na(idf$SMWT)),]                             
  
  # Subset the baseline data as training set                                    ## 231 cases remain
  train <- idf[idf$Visit == "V2",]
  train$Visit <- NULL
  
  # Box-Cox-transformation of dependent variables
  train$SMWT <- boxcox_transformation(train$SMWT)
  train$Mode <- boxcox_transformation(train$Mode)
  
  # Fit cv.glmnet function on bootstrapped training datasets                    
  boot_res <- boot(data=train, fitglmnet, R=5000, alpha=0.5, 
                   nfolds=10, dvariable=dvariables, family="mgaussian")
  
  # Obtain results and calculate variable inclusion frequencies (VIP's)
  boot_coef <- as.data.frame(boot_res$t)
  colnames(boot_coef) <- names(boot_res$t0)                            
  boot_freq <- round((colSums(boot_coef != 0) / 5000),2)
  
  # Add results to imputation overarching dataframe
  selected_pred <- bind_rows(selected_pred, boot_freq)
  
}

## Calculate average VIP and add to VIP dataframe         
cm <- colMeans(selected_pred)
selected_pred <- rbind(selected_pred, cm)
row.names(selected_pred) <- c(1:10, "Mean")

## Transpose, order and restore protein ID's
t_selected_pred <- as.data.frame(t(selected_pred))
t_selected_pred <- t_selected_pred[order(t_selected_pred$Mean, decreasing = T),]

for (i in 1:nrow(t_selected_pred)){
  prot <- rownames(t_selected_pred)[i]
  rownames(t_selected_pred)[i] <- 
    namekey$Original[namekey$Permuted==prot]
}
  
## Save results
save(cdf_e, cdf, t_selected_pred, namekey,
     file = "p_DIA_NN_CTG&SMWT_prot_VIPs.RDATA")
write.csv(t_selected_pred, file="CTG&SMWT_MVpred_VIPs.csv")


#####################################################################
## Internal validation of multivariate predictors for CTG and SMWT ##
#####################################################################

## Optional: Load results (manually reload libraries and box-cox function from above)
load(file = "p_DIA_NN_CTG&SMWT_prot_VIPs.RDATA")

## Generate lists to store linear models
rec_MV_fits <- list()
rec_CTG_fits <- list()
rec_SMWT_fits <- list()

## Generate dataframes to store statistical performance measures
CTG_stat_perf <- data.frame(matrix(ncol=9, nrow=0))
colnames(CTG_stat_perf) <- c("VIP (%)", "# predictors",
                             "Train_RSQ", "Train_aRSQ", "Train_RMSE",
                             "Test_RSQ", "Test_RMSE",
                             "Val_RSQ", "Val_RMSE")
SMWT_stat_perf <- CTG_stat_perf

## Change protein-group names according to namekey for linear model fits
# Note: The lm() function does not accept a protein group name as predictor
for (c in 1:nrow(t_selected_pred)){
  rownames(t_selected_pred)[c] <- 
    namekey$Permuted[namekey$Original==rownames(t_selected_pred)[c]]
}

for (PID in namekey$Original){
  colnames(cdf)[colnames(cdf)==PID] <-
    namekey$Permuted[namekey$Original==PID]
}

for (PID in namekey$Original){
  if (PID %in% colnames(cdf_e)){
     colnames(cdf_e)[colnames(cdf_e)==PID] <-
       namekey$Permuted[namekey$Original==PID]
  }
}

## Generate results for different variable inclusion probabilities (VIPs)
for (i in 9:5){
  i <- i/10
  a <- nrow(CTG_stat_perf) + 1
  
  ## Obtain predictors with VIP according to i and amount to dataframe
  CTG_SMWT_predictors <- rownames(t_selected_pred[t_selected_pred$Mean >= i,])
  CTG_stat_perf[a, "VIP (%)"] <- i*100
  CTG_stat_perf[a, "# predictors"] <- length(CTG_SMWT_predictors)
  SMWT_stat_perf[a, "VIP (%)"] <- i*100
  SMWT_stat_perf[a, "# predictors"] <- length(CTG_SMWT_predictors)
  
  ## Subset relevant data, split training & testing
  cdf_pred <- cdf[,colnames(cdf) %in% c("SMWT", "Mode", "Visit", CTG_SMWT_predictors)]
  
  # Split training & testing set
  train <- cdf_pred[cdf_pred$Visit == "V2",]
  test <- cdf_pred[cdf_pred$Visit == "V4",] 
  
  ## Transform and Standardize dependent variables
  train$Mode <- scale(boxcox_transformation(train$Mode))[1:nrow(train)]
  train$SMWT <- scale(boxcox_transformation(train$SMWT))[1:nrow(train)]
  test$Mode <- scale(boxcox_transformation(test$Mode))[1:nrow(test)]
  test$SMWT <- scale(boxcox_transformation(test$SMWT))[1:nrow(test)]
  
  ## Standardize independent variables
  train[,4:ncol(train)] <- scale(train[,4:ncol(train)])
  test[,4:ncol(test)] <- scale(test[,4:ncol(test)])
  
  ## Build Multivariate linear model & add to list
  fitname <- paste("VIP%", i)
  formula = paste("cbind(Mode, SMWT)", paste(CTG_SMWT_predictors, collapse = " + "), sep=" ~ ")
  mv_fit <- lm(formula=formula, data=train)
  rec_MV_fits[[fitname]] <- mv_fit
  
  ## Build individual linear models & add to list [for model assumption validations]
  CTG_formula = paste("Mode", paste(CTG_SMWT_predictors, collapse = " + "), sep=" ~ ")
  rec_CTG_fits[[fitname]] <- lm(formula=CTG_formula, data=train)
  SMWT_formula = paste("SMWT", paste(CTG_SMWT_predictors, collapse = " + "), sep=" ~ ")
  rec_SMWT_fits[[fitname]] <- lm(formula=SMWT_formula, data=train)
  
  ## Predict training and testing data & add to respective dataframes
  train_pred <- as.data.frame(predict(mv_fit, newdata=train))
  colnames(train_pred) <- c("CTG_pred", "SMWT_pred")
  test_pred <- as.data.frame(predict(mv_fit, newdata=test))
  colnames(test_pred) <- c("CTG_pred", "SMWT_pred")
  
  table(rownames(train_pred) == rownames(train))
  train <- merge(train, train_pred, by="row.names")
  table(rownames(test_pred) == rownames(test))
  test <- merge(test, test_pred, by="row.names")
  
  ## Predict external validation data
  # Subset relevant data in predictor dataframe & Transform and standardize dependent variables
  cdf_e_pred <- cdf_e[,colnames(cdf_e) %in% c("CTG", "SMWT", CTG_SMWT_predictors)]
  cdf_e_pred$CTG <- scale(boxcox_transformation(cdf_e_pred$CTG))[1:nrow(cdf_e_pred)]
  cdf_e_pred$SMWT <- scale(boxcox_transformation(cdf_e_pred$SMWT))[1:nrow(cdf_e_pred)]
  
  # Standardize independent variables
  cdf_e_pred[,3:ncol(cdf_e_pred)] <- scale(cdf_e_pred[,3:ncol(cdf_e_pred)])
  
  # Predict CTG and SMWT scores based on ReCognitION training data model
  cdf_e_pred$CTG_predicted <- predict(mv_fit, newdata=cdf_e_pred)[,1]
  cdf_e_pred$SMWT_predicted <- predict(mv_fit, newdata=cdf_e_pred)[,2]
  
  ## Obtain / Calculate performance metrics for CTG and add to dataframe
  # Training data
  CTG_stat_perf[a, "Train_RSQ"] <- round(summary(mv_fit)$'Response Mode'$r.squared, 4)
  CTG_stat_perf[a, "Train_aRSQ"] <- round(summary(mv_fit)$'Response Mode'$adj.r.squared, 4)
  CTG_train <- train[,c("Mode", "CTG_pred")]
  CTG_train <- CTG_train[complete.cases(CTG_train),]
  CTG_stat_perf[a, "Train_RMSE"] <- round(rmse(CTG_train$Mode, CTG_train$CTG_pred), 4)
  
  # Testing data
  CTG_test <- test[,c("Mode", "CTG_pred")]
  CTG_test <- CTG_test[complete.cases(CTG_test),]
  CTG_stat_perf[a, "Test_RSQ"] <- round(rsq(CTG_test$Mode, CTG_test$CTG_pred), 4)
  CTG_stat_perf[a, "Test_RMSE"] <- round(rmse(CTG_test$Mode, CTG_test$CTG_pred), 4)
  
  # Validation data
  cdf_e_CTG <- cdf_e_pred[complete.cases(cdf_e_pred[,c("CTG", "CTG_predicted")]),]
  CTG_stat_perf[a, "Val_RSQ"] <- round(rsq(cdf_e_CTG$CTG, cdf_e_CTG$CTG_predicted),4)
  CTG_stat_perf[a, "Val_RMSE"] <- round(rmse(cdf_e_CTG$CTG, cdf_e_CTG$CTG_predicted),4)
  
  ## Obtain / Calculate performance metrics for SMWT and add to dataframe
  # Training data
  SMWT_stat_perf[a, "Train_RSQ"] <- round(summary(mv_fit)$'Response SMWT'$r.squared, 4)
  SMWT_stat_perf[a, "Train_aRSQ"] <- round(summary(mv_fit)$'Response SMWT'$adj.r.squared, 4)
  SMWT_train <- train[,c("SMWT", "SMWT_pred")]
  SMWT_train <- SMWT_train[complete.cases(SMWT_train),]
  SMWT_stat_perf[a, "Train_RMSE"] <- round(rmse(SMWT_train$SMWT, SMWT_train$SMWT_pred),4)
  
  # Testing data
  SMWT_test <- test[,c("SMWT", "SMWT_pred")]
  SMWT_test <- SMWT_test[complete.cases(SMWT_test),]
  SMWT_stat_perf[a, "Test_RSQ"] <- round(rsq(SMWT_test$SMWT, SMWT_test$SMWT_pred),4)
  SMWT_stat_perf[a, "Test_RMSE"] <- round(rmse(SMWT_test$SMWT, SMWT_test$SMWT_pred),4)
  
  # Validation data
  cdf_e_SMWT <- cdf_e_pred[complete.cases(cdf_e_pred[,c("SMWT", "SMWT_predicted")]),]
  SMWT_stat_perf[a, "Val_RSQ"] <- round(rsq(cdf_e_SMWT$SMWT, cdf_e_SMWT$SMWT_predicted),4)
  SMWT_stat_perf[a, "Val_RMSE"] <- round(rmse(cdf_e_SMWT$SMWT, cdf_e_SMWT$SMWT_predicted),4)
}

## Combine CTG and SMWT results for different VIP's and save
CTG_stat_perf$Outcome <- "CTG"
CTG_stat_perf <- CTG_stat_perf[, c(10, 1:9)]
SMWT_stat_perf$Outcome <- "SMWT"
SMWT_stat_perf <- SMWT_stat_perf[,c(10,1:9)]
table(colnames(CTG_stat_perf)==colnames(SMWT_stat_perf))
T4 <- rbind(CTG_stat_perf, SMWT_stat_perf)
write.csv(T4, file="T4_VIP_performance_comp.csv")

## Obtain optimal models
rec_MV_fit <- rec_MV_fits[[4]]
rec_CTG_fit <- rec_CTG_fits[[4]]
rec_SMWT_fit <- rec_SMWT_fits[[4]]

## Check assumptions for individual univariate linear regressions 
gvlma(rec_CTG_fit)   # No pass on Link Function.
gvlma(rec_SMWT_fit)  # All assumptions acceptable.

# Run diagnostics on CTG fit model
plot(rec_CTG_fit, which=4) #Outlier A022_V2

## Remove outlier from datasets & refit models
MV_formula <- formula(rec_MV_fit)
MV_depvar <- as.data.frame(model.response(rec_MV_fit$model))
MV_train <- cbind(MV_depvar, rec_MV_fit$model[,2:14])
MV_train <- MV_train[rownames(MV_train) != "A022_V2",]
rec_MV_fit <- lm(formula=MV_formula, data=MV_train)

CTG_formula <- formula(rec_CTG_fit)
CTG_train <- rec_CTG_fit$model[rownames(rec_CTG_fit$model) != "A022_V2",]
rec_CTG_fit <- lm(formula=CTG_formula, data=CTG_train)

SMWT_formula <- formula(rec_SMWT_fit)
SMWT_train <- rec_SMWT_fit$model[rownames(rec_SMWT_fit$model) != "A022_V2",]
rec_SMWT_fit <- lm(formula=SMWT_formula, data=SMWT_train)

## Re-check assumptions for individual univariate linear regressions 
gvlma(rec_CTG_fit)   # All assumptions acceptable.
gvlma(rec_SMWT_fit)  # All assumptions acceptable. 

## Generate STable 1: ReCognitION multivariate fit coefficients and p-values & save
rec_MV_fit_sum <- summary(rec_MV_fit)
T5_CTG <- as.data.frame(rec_MV_fit_sum$`Response Mode`$coefficients)[,c(1,2,4)]
colnames(T5_CTG) <- c("CTG_estimate", "CTG_Std.Error", "CTG_pvalue")

T5_SMWT <- as.data.frame(rec_MV_fit_sum$`Response SMWT`$coefficients)[,c(1,2,4)]
colnames(T5_SMWT) <- c("SMWT_estimate", "SMWT_Std.Error", "SMWT_pvalue")

T5 <- as.data.frame(cbind(T5_CTG, T5_SMWT))
T5 <- apply(T5, 2, function(x){
  round(x, digits = 4)
})

for (i in 2:nrow(T5)){
  rownames(T5)[i] <- namekey$Original[namekey$Permuted==rownames(T5)[i]]
}

write.csv(T5, file="T5_MV_CTG_SMWT_pred_rec_coef.csv")
