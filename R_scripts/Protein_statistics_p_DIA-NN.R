
## Script generated on 26/1/23 by D. van As
## Purpose: Statistical analysis of pooled DIA-NN Protein data
## Last update: 10/06/24


####################
## Load libraries ##
####################
library("lme4")     #V 1.1.29
library("lmerTest") #V 3.1.3


################
## Load files ##
################

## Filtered and normalized protein counts, phenotype data
load(file = "normalized_protein_intensities_p_DIA_NN.RDATA")

# Add patient ID's
pdf$PatientID <- gsub("_V.", "", rownames(pdf))

## Combine protein counts and phenotype data
rpc_t <- as.data.frame(t(rpc))
table(rownames(rpc_t)==rownames(pdf))
cdf <- merge(pdf, rpc_t, by="row.names")
rownames(cdf) <- cdf$Row.names
cdf$Row.names <- NULL


###################################################################### 
## Model 2: CTG-repeat association, corrected for Visit, Sex, Plate ##
######################################################################

## Control: Randomize CTG repeat
#set.seed(1)
#cdf$Mode <- sample(cdf$Mode)

## Fit Mixed effect models
fit_list <- list()
for (protein in colnames(cdf)[37:ncol(cdf)]){                                   
  fit <- suppressMessages(
    lmer(formula = cdf[,protein] ~ SexCode + Visit + Mode + (1|plate) + (1|PatientID), data=cdf))
  ## Exclude fits with convergence problems
  sfit <- summary(fit)
  if (is.null(sfit$optinfo$conv$lme4$code)){
    fit_list[[protein]] <- fit
  }
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
}
Mode_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][4,]})))

## Add FDR correction & specify column names
Mode_res_df$FDR <- p.adjust(Mode_res_df$`Pr(>|t|)`, method='fdr')
colnames(Mode_res_df) <- paste("Mode_", colnames(Mode_res_df), sep="" )
rownames(Mode_res_df) <- names(fit_coefficients)

## Add correlation coefficients
cor_vec <- c()
for (protein in rownames(Mode_res_df)){
  cor_vec <- c(cor_vec, cor(cdf[,protein], cdf[,"Mode"], use="pairwise.complete.obs", method="pearson"))
}
Mode_res_df$pcor <- cor_vec

## Check for significant hits                                                   
sum(Mode_res_df$Mode_FDR < 0.05)                                                # 161 hits; Randomized CTG: 0 hits


###########################################################################                                     
## Models 3: Individual outcome measures, Sex, Visit and Plate corrected ##
###########################################################################

## Control: Randomize SMWT                                                      # only 1 FDR hit after randomization       
#set.seed(1)
#cdf$SMWT <- sample(cdf$SMWT)

## Run analysis for each outcome measure, store results in list outcome_coef
outcome_coef <- list()

for (OM in colnames(cdf)[c(1:2, 4:10, 12:20, 22, 25:30)]){                                                
  
  ## Print status
  print(paste("Currently fitting ", OM))
  
  ## Fit Mixed effect models
  fit_list <- list()
  for (protein in colnames(cdf)[37:ncol(cdf)]){                                                                 
    ## Model fit
    fit <- suppressMessages(
      lmer(formula = cdf[,protein] ~ SexCode + Visit + cdf[,OM] + (1|plate) + (1|PatientID), data=cdf))
    ## Exclude fits with convergence problems
    sfit <- summary(fit)
    if (is.null(sfit$optinfo$conv$lme4$code)){
      fit_list[[protein]] <- fit
    }
  }
  
  ## Obtain coefficient estimates and p-values
  fit_coefficients <- list()
  for (fit in names(fit_list)){
    fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
  }
  
  ## Split and merge results per predictor
  OM_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][4,]})))
  
  ## Add FDR correction & specify column names
  OM_res_df$FDR <- p.adjust(OM_res_df$`Pr(>|t|)`, method='fdr')
  colnames(OM_res_df) <- paste(OM, colnames(OM_res_df), sep="_" )
  rownames(OM_res_df) <- names(fit_coefficients)
  
  ## Add correlation coefficient if OM is numeric
  if (class(cdf[,OM]) == "numeric"){
    cor_vec <- c()
    for (protein in rownames(OM_res_df)){
      cor_vec <- c(cor_vec, cor(cdf[,protein], cdf[,OM], use="pairwise.complete.obs", method="pearson"))
    }
    OM_res_df$pcor <- cor_vec
  }
  
  ## Store results
  outcome_coef[[OM]] <- OM_res_df
}


## Overview of significant associations per outcome measure                     
coef_table = data.frame()

for (i in 1:length(outcome_coef)){
  # Count number of significant hits per predictor
  num_p <- length(outcome_coef[[i]][,3][outcome_coef[[i]][,3] < 0.05]) 
  num_FDR <- length(outcome_coef[[i]][,4][outcome_coef[[i]][,4] < 0.05])
  
  # Store in table
  num_df <- cbind(num_p, num_FDR)
  coef_table <- rbind(coef_table, num_df)
}
rownames(coef_table) <- names(outcome_coef)
coef_table <- coef_table[order(coef_table$num_FDR, decreasing = T),]


###############################################################                    
## Model 4: CBT Intervention effect Sex and Plate corrected ##
###############################################################

## Fit Mixed effect models
fit_list <- list()
for (protein in colnames(cdf)[37:ncol(cdf)]){                                   
  fit <- suppressMessages(
    lmer(formula = cdf[,protein] ~ SexCode + Visit*TreatmentCode + (1|plate) + (1|PatientID), data=cdf))
  ## Exclude fits with convergence problems
  sfit <- summary(fit)
  if (is.null(sfit$optinfo$conv$lme4$code)){
    fit_list[[protein]] <- fit
  }
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
}
Vis_Treat_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][5,]})))

## Add FDR correction & specify column names
Vis_Treat_res_df$FDR <- p.adjust(Vis_Treat_res_df$`Pr(>|t|)`, method='fdr')
colnames(Vis_Treat_res_df) <- paste("Visit_Treatment_", colnames(Vis_Treat_res_df), sep="" )
rownames(Vis_Treat_res_df) <- names(fit_list)

## Check for significant hits
sum(Vis_Treat_res_df$Visit_Treatment_FDR < 0.05)
Vis_Treat_res_df[Vis_Treat_res_df$Visit_Treatment_FDR < 0.05, ]

## Plot p-value distributions of predictors
hist(Vis_Treat_res_df$'Visit_Treatment_Pr(>|t|)', n=50,
     xlab="P-values",
     main="Protein: P-value distribution of Treatment*Visit effect")

#### Supplemental analysis: no plate correction

## Fit Mixed effect models
fit_list <- list()
for (protein in colnames(cdf)[37:ncol(cdf)]){                                   
  fit <- suppressMessages(
    lmer(formula = cdf[,protein] ~ SexCode + Visit*TreatmentCode + (1|PatientID), data=cdf))
  ## Exclude fits with convergence problems
  sfit <- summary(fit)
  if (is.null(sfit$optinfo$conv$lme4$code)){
    fit_list[[protein]] <- fit
  }
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
}
Vis_treat_sup_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][5,]})))

## Add FDR correction & specify column names
Vis_treat_sup_df$FDR <- p.adjust(Vis_treat_sup_df$`Pr(>|t|)`, method='fdr')
colnames(Vis_treat_sup_df) <- paste("Visit_Treatment_", colnames(Vis_treat_sup_df), sep="" )
rownames(Vis_treat_sup_df) <- names(fit_list)

## Comparison of results
table(rownames(Vis_treat_sup_df[Vis_treat_sup_df$Visit_Treatment_FDR < 0.05,]) %in% 
  rownames(Vis_Treat_res_df[Vis_Treat_res_df$Visit_Treatment_FDR < 0.05,]))
# Conclusion supplemental analysis: no new hits without plate effect correction. 


############################################################################                     
## Model 5: Graded Exercise Therapy (GET) effect Sex and Plate corrected ##
############################################################################

## Fit Mixed effect models
fit_list <- list()
for (protein in colnames(cdf)[37:ncol(cdf)]){                                   
  fit <- suppressMessages(
    lmer(formula = cdf[,protein] ~ SexCode + Visit*GradedExerciseTherapy + (1|plate) + (1|PatientID), data=cdf))
  ## Exclude fits with convergence problems
  sfit <- summary(fit)
  if (is.null(sfit$optinfo$conv$lme4$code)){
    fit_list[[protein]] <- fit
  }
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
}
Vis_Get_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][5,]})))

## Add FDR correction & specify column names
Vis_Get_res_df$FDR <- p.adjust(Vis_Get_res_df$`Pr(>|t|)`, method='fdr')
colnames(Vis_Get_res_df) <- paste("Visit_Get_", colnames(Vis_Get_res_df), sep="" )
rownames(Vis_Get_res_df) <- names(fit_list)

## Check for significant hits
sum(Vis_Get_res_df$Visit_Get_FDR < 0.05)
Vis_Get_res_df <- Vis_Get_res_df[order(Vis_Get_res_df$`Visit_Get_Pr(>|t|)`),]

## Plot p-value distributions of predictors
hist(Vis_Get_res_df$`Visit_Get_Pr(>|t|)`, n=50,
     xlab="P-values",
     main="Protein: P-value distribution of Treatment*Get effect")

#### Supplemental analysis: no plate correction

## Fit Mixed effect models
fit_list <- list()
for (protein in colnames(cdf)[37:ncol(cdf)]){                                   
  fit <- suppressMessages(
    lmer(formula = cdf[,protein] ~ SexCode + Visit*GradedExerciseTherapy + (1|PatientID), data=cdf))
  ## Exclude fits with convergence problems
  sfit <- summary(fit)
  if (is.null(sfit$optinfo$conv$lme4$code)){
    fit_list[[protein]] <- fit
  }
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,5)]
}
Vis_Get_sup_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][5,]})))

## Add FDR correction & specify column names
Vis_Get_sup_df$FDR <- p.adjust(Vis_Get_sup_df$`Pr(>|t|)`, method='fdr')
colnames(Vis_Get_sup_df) <- paste("Visit_Get_", colnames(Vis_Get_sup_df), sep="" )
rownames(Vis_Get_sup_df) <- names(fit_list)

## Comparison of results
table(rownames(Vis_Get_sup_df[Vis_Get_sup_df$Visit_Get_FDR < 0.05,]) %in%
  rownames(Vis_Get_res_df[Vis_Get_res_df$Visit_Get_FDR < 0.05,]))


##################
## Save results ##
##################

## Save results
save(cdf, Vis_Treat_res_df, Vis_Get_res_df, Mode_res_df, outcome_coef, coef_table,
     file = "DIA-NN_p_Protein_statistics.RDATA")

## Generate excel output 
SMWT_df <- outcome_coef[["SMWT"]]

write.csv(Mode_res_df, file="ReCognitION_CTG_hits.csv")
write.csv(SMWT_df, file="ReCognitION_SMWT_hits.csv")
write.csv(Vis_Treat_res_df, file="ReCognitION_CBT_hits.csv")
write.csv(Vis_treat_sup_df, file="ReCognitION_CBT_no_plate_cor_hits.csv")
write.csv(Vis_Get_res_df, file="ReCognitION_GET_hits.csv")
write.csv(Vis_Get_sup_df, file="ReCognitION_GET_no_plate_cor_hits.csv")
write.csv(coef_table, file="ReCognitION_OM_summary.csv")
