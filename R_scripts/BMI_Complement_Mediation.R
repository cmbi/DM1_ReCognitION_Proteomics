
## Script generated on 13.06.25, last update on: 18/11/25
## Mediation analysis to see if complement activation is moderated by BMI

###############
## Libraries ##
###############
library(readxl)    #V 1.4.5
library(mediation) #V 4.5.1


############################
## Load relevant datasets ##
############################

## ReCognitION proteomics statistical results
load (file = "DIA-NN_p_Protein_statistics.RDATA")

## Table 2 to obtain overlapping significant/validated protein groups
load(file = "Table2.RDATA")

## BMI data from OPTIMISTIC
# Post-hoc, peer-review requested analysis, BMI dropped in original phenotype data
BMI_df <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=1)
BMI_df <- as.data.frame(BMI_df)


##################################
## Subset & merge relevant data ##
##################################

# Height & Weight data only recorded at trial start, analyses therefore limited to V2

## Subset V2 data of proteomics dataframe
cdf_v2 <- cdf[grepl("_V2", rownames(cdf)),]

## Subset BMI data
BMI_df <- BMI_df[,c("PatientID", "Height", "Weight")]
BMI_df$PatientID <- gsub("P", "", BMI_df$PatientID)
BMI_df <- BMI_df[BMI_df$PatientID %in% cdf_v2$PatientID,]

## Merge
table(cdf_v2$PatientID %in% BMI_df$PatientID)
table(BMI_df$PatientID %in% cdf_v2$PatientID)
nrow(cdf_v2) == nrow(BMI_df)
df <- merge(cdf_v2, BMI_df, by="PatientID")


####################
## Pre-processing ##
####################

## Check for outliers 
sort(df$Height)
df$Height[df$Height > 10] <- NA                                                 ## Two cases values removed
sort(df$Weight)

## Calculate BMI
df$BMI <- df$Weight/(df$Height^2)
hist(df$BMI, n=30, xlab="BMI", main="OPTIMISTIC baseline cohort")               ## Two very high BMI cases, not removed


########################
## Mediation analysis ##
########################

options(scipen = 999)

## Set up dataframe to store results
ST1 <- data.frame()

## CTG-repeat mediation analysis
# Define variables and create vector to store results
X <- "Mode"          # Covariate of interest (Mode = CTG-repeat)
M <- "BMI"           # Mediator (BMI)
Y <- ""              # Protein expression

for (Y in T2$Rec_PG[T2$OM=="CTG"]){
  
  # Subset relevant data and select complete cases
  df_CTG <- df[, c(Y, X, M)]
  df_CTG <- df_CTG[complete.cases(df_CTG),]

  # Fit mediator model: How CTG-repeat affects BMI
  formula.M <- paste(M, X, sep=" ~ ")
  model.M <- lm(formula.M, data=df_CTG)
  
  # Fit outcome model: How CTG-repeat and BMI affect protein expression
  formula.Y <- paste("`", Y, "`", " ~ ", X, " + ", M, sep='')
  model.Y <- lm(formula.Y, data=df_CTG)
  
  # Run mediation analysis and obtain results
  MRes <- mediate(model.M, model.Y, treat = X, mediator = M,                    
                  boot=T, sims=1000)
  ACME <- round(MRes$d0, 6)
  ADE <- round(MRes$z0, 6)
  TOT <- round(MRes$tau.coef, 6)
  prop_med <- round(MRes$d0 / MRes$tau.coef, 6)
  
  # Store all results in dataframe
  res <- c(Y, ACME, ADE, TOT, prop_med)
  ST1 <- rbind(ST1, res)
}

## 6MWT mediation analysis
# Define variables and create vector to store results
X <- "SMWT"           # Covariate of interest (6MWT)
M <- "BMI"            # Mediator (BMI)
Y <- ""               # Protein expression

for (Y in T2$Rec_PG[T2$OM=="6MWT"]){
  
  # Subset relevant data and select complete cases
  df_smwt <- df[, c(Y, X, M)]
  df_smwt <- df_smwt[complete.cases(df_smwt),]
  
  # Fit mediator model: How SMWT affects BMI
  formula.M <- paste(M, X, sep=" ~ ")
  model.M <- lm(formula.M, data=df_smwt)
  
  # Fit outcome model: How SMWT and BMI affect protein expression
  formula.Y <- paste("`", Y, "`", " ~ ", X, " + ", M, sep='')
  model.Y <- lm(formula.Y, data=df_smwt)
  
  # Run mediation analysis and obtain results
  MRes <- mediate(model.M, model.Y, treat = X, mediator = M,                  
                  boot=T, sims=1000)
  ACME <- round(MRes$d0, 6)
  ADE <- round(MRes$z0, 6)
  TOT <- round(MRes$tau.coef, 6)
  prop_med <- round(MRes$d0 / MRes$tau.coef, 6)
  
  # Store all results in dataframe
  res <- c(Y, ACME, ADE, TOT, prop_med)
  ST1 <- rbind(ST1, res)
}

## Correct column names and add co-variate information
predictor <- c(rep("CTG-repeat", length(T2$Rec_PG[T2$OM=="CTG"])),
           rep("6MWT", length(T2$Rec_PG[T2$OM=="6MWT"])))
ST1 <- cbind(ST1, predictor)

colnames(ST1) <- c("Protein(group)", "ACME", "ADE", "Total", "Prop. Mediated", "Predictor")
ST1 <- ST1[,c(6, 1, 2, 3, 4, 5)]


## Save results
write.csv(ST1, file="STable1.csv")
save(ST1, 
     file = "STable1.RDATA")








