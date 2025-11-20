

## Script generated on 19/8/22 by D. van As
## Last updates applied on: 18/11/25: phenotype data contains multiple tabs now; updated for easy reproduciblity 
## Purpose: Calculate combined scores, filter out certain non-continuous outcomes, remove outliers


###############
## Libraries ##
###############

library(readxl) #V 1.4.0


###############                                                                 
## Functions ##
###############

IQR_filter <- function(x, n){
  # This function removes outliers from a vector based on the Interquartile Range Rule
  # The upper and lower boundaries are respectively Q1-n*IQR and Q3+n*IQR
  Q <- quantile(x, probs=c(0.25, 0.75), na.rm = T)
  iqr <- IQR(x, na.rm = T)
  low <- unname(Q[1] - n*iqr)
  up <- unname(Q[2] + n*iqr)
  x[x < low | x > up] <- NA
  return(x)
}


###########################
## Load relevant datasets##
###########################

## Clinical trial data 
dfa <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=1)
dfb <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=2)
dfc <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=3)
dfd <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=4)
dfe <- read_excel("OPTIMISTIC_Phenotype_Data.xlsx", sheet=5)
pdf <- cbind(dfa, dfb, dfc, dfd, dfe)
rownames(pdf) <- pdf$PatientID


######################################
## Calculate stroop, TMT and V4Age  ##
######################################

# Stroop Interference 
stroop2V2 <- ((60-pdf$StroopCardIIErrorsV2) / 60) / pdf$StroopCardIITimeV2
stroop2V4 <- ((60-pdf$StroopCardIIErrorsV4) / 60) / pdf$StroopCardIITimeV4
stroop3V2 <- ((60-pdf$StroopCardIIIErrorsV2) / 60) / pdf$StroopCardIIITimeV2
stroop3V4 <- ((60-pdf$StroopCardIIIErrorsV4) / 60) / pdf$StroopCardIIITimeV4
pdf$StroopInterferenceV2 <- stroop3V2 / stroop2V2
pdf$StroopInterferenceV4 <- stroop3V4 / stroop2V4

# Trail Making Test
pdf$TMTV2 <- pdf$TMTBV2 / pdf$TMTAV2
pdf$TMTV4 <- pdf$TMTBV4 / pdf$TMTAV4

## Check & correct graded-exercise therapy versus intervention group
# Patients of the control group did not receive graded exercise therapy
table(pdf$GradedExerciseTherapy[pdf$TreatmentCode==0])
pdf$PatientID[pdf$TreatmentCode==0 & pdf$GradedExerciseTherapy==1]              #C024P
pdf$GradedExerciseTherapy[pdf$PatientID == "C024P"] <- 0


######################################################
## Subset relevant V2 and V4 measurements and merge ##
######################################################

## V2 (= Baseline) Measurements
pdfV2 <- pdf[,grep("V2", colnames(pdf))]
pdfV2$Visit <- "V2"
pdfV2 <- cbind(pdf[,c(3, 5, 8, 10, 27)], pdfV2)
pdfV2 <- pdfV2[,-c(7:8, 10, 19:20, 24:29, 34)]
colnames(pdfV2) <- gsub("V2", "", colnames(pdfV2))
rownames(pdfV2) <- paste(rownames(pdfV2), pdfV2$Visit, sep="")

## V4 (= 10 months) Measurements
pdfV4 <- pdf[,grep("V4", colnames(pdf))]
pdfV4$Mode <- pdf$V2Mode
pdfV4$Visit <- "V4"
pdfV4$Age <- pdf$V2Age+0.83                                                     
pdfV4 <- cbind(pdf[,c(3, 5, 8, 10, 27)], pdfV4)
pdfV4 <- pdfV4[,-c(14:15, 19:24)]
colnames(pdfV4) <- gsub("V4", "", colnames(pdfV4))
rownames(pdfV4) <- paste(rownames(pdfV4), pdfV4$Visit, sep="")

## Merge & correct rownames
pdfV2 <- pdfV2[,order(colnames(pdfV2))]
pdfV4 <- pdfV4[,order(colnames(pdfV4))]
table(colnames(pdfV2) == colnames(pdfV4))
pdf <- rbind(pdfV2, pdfV4)                                                      
rownames(pdf) <- gsub("P", "_", rownames(pdf))

## Drop non-relevant or non-continuous outcome measures
# Variables dropped: Borg, MIRS
pdf <- pdf[,-c(21,23)]


####################################
## Correct variable types & names ##
####################################                                            

# Correct variable types                                                        
factor <- c("TreatmentCode", "GradedExerciseTherapy", "SexCode",
            "SiteCode", "VariantRepeats", "Visit")
numeric <- colnames(pdf)[!colnames(pdf)%in% c(factor)]

for (variable in colnames(pdf)){
  if (variable %in% numeric){
    pdf[,variable] <- as.numeric(as.character(pdf[,variable]))
  }
  if (variable %in% factor){
    pdf[,variable] <- as.factor(as.character(pdf[,variable]))
  }
}

# Confirm variable types
lapply(pdf, class)

# Correct some variable names
colnames(pdf)[colnames(pdf)=="INQOLQolScore"] <- "INQoL"
numeric[13] <- "INQoL"

######################
## Outlier handling ##
######################

cdf <- pdf #Generate backup

## Check NA's across whole dataset
sum(is.na(cdf))                                                                 # 2301 NA 

## Apply IQR filter to numeric/continuous variables
for(OM in numeric){                                                             # 27 numeric variables *  510 rows = 13770
  cdf[,OM] <- IQR_filter(cdf[,OM], 3)
}
sum(is.na(cdf))                                                                 # (10 values flagged / (13770-2301))*100 = 0.09%                                                        

## Manually check flagged outliers
apply(pdf, 2, function(x) length(which(!is.na(x)))) - apply(cdf, 2, function(x) length(which(!is.na(x))))

## Accept all changes
pdf <- cdf


##################################
## Save filtered phenotype data ##
##################################

save(pdf, 
     file= "Filtered_Phenotypedata.RData")


