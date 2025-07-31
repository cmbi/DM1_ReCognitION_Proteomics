
## Script generated on 23/1/23 by D van As
## Purpose of script: ReCognitION DIA-pooled DIA-NN protein sample pre-filter analysis
## Last updates on: 02/06/25


####################
## Load libraries ##
####################
library("readxl")        #V 1.4.0
library("DEqMS")         #V 1.8.0
library("RColorBrewer")  #V 1.1.3
library("mice")          #V 3.14.0
library("micemd")        #V 1.8.0


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


############################
## Load relevant datasets ##
############################

## Load ID_vec from peptide results
load(file = "normalized_peptide_intensities_p_DIA_NN.RData")
rm(pdf, rpc)

## Raw Protein counts
rpc <- read.table("merged_speclib.pg_matrix.tsv", sep = "\t", header=T)

## Sample randomization key
pkey <- read_excel("All samples_LOG_randomize_runorder.xlsx", sheet=2)
pkey <- as.data.frame(pkey)
pkey <- pkey[-c(464:468),] # Reserve samples, pools

## Phenotype data
# Results generated using Phenotype_filter_imputation.R
load(file= "Filtered_Phenotypedata.RDATA")


#########################
## Data pre-processing ##
#########################

## Remove spike-in bovine-albumin and pig-trypsin
any(grepl("P00764", rpc$Protein.Group))
any(grepl("TRYP_PIG", rpc$Protein.Group))
rpc[grepl("TRYP_PIG", rpc$Protein.Group),1:5]
rpc <- rpc[!grepl("TRYP_PIG", rpc$Protein.Group),]

any(grepl("P02769", rpc$Protein.Group))
any(grepl("ALBU_BOVIN", rpc$Protein.Group))
rpc[grepl("ALBU_BOVIN", rpc$Protein.Group),1:5]
rpc <- rpc[!grepl("ALBU_BOVIN", rpc$Protein.Group),]

## Simplify protein sample names to randomization numbers
colnames(rpc)[6:ncol(rpc)] <- gsub(".*iRT_", "", colnames(rpc)[6:ncol(rpc)])
colnames(rpc)[colnames(rpc) == "2G6_20210817142613.raw"] <- "2G6.raw"
colnames(rpc)[colnames(rpc) == "1D1_20210820165914.raw"] <- "1D1.raw"
colnames(rpc)[colnames(rpc) == "3F12_20210830074845.raw"] <- "3F12.raw"
colnames(rpc)[6:ncol(rpc)] <- gsub(".raw", "", colnames(rpc)[6:ncol(rpc)])

## Change 0's to NA's 
rpc[rpc == 0] <- NA

## Set protein group names as row names & drop protein, NumPeptides and PeptideSequences columns
rownames(rpc) <- rpc$Protein.Group
rpc <- rpc[,-c(1:5)]

## Remove samples based on peptide analysis                                     # 10 samples removed due to technical problems / low total counts
rpc[,!colnames(rpc) %in% ID_vec] <- NULL

## Set intensities of < 10.000 to NA (unreliable)
rpc[rpc < 10000] <- NA


#########################################                    
## Log2 transformation of protein data ##
#########################################

## Check if all rpcs are numeric | All intensities stored as numeric
table(apply(rpc, 2, function(x) is.numeric(x)))

## Log2(x) transformation of intensities
rpc <- log2(rpc)


########################################                                        
## Median based sample normalization  ##
########################################

# Box plots before scaling
boxplot(rpc[,1:150], ylab="Protein intensities [log2]")
boxplot(rpc[,151:300], ylab="Protein intensities [log2]")
boxplot(rpc[,301:ncol(rpc)], ylab="Protein intensities [log2]")

# Equal Median Normalization
rpc <- equalMedianNormalization(rpc) 

# Box plots after scaling
boxplot(rpc[,1:150], ylab="Protein intensities [log2]")                         
boxplot(rpc[,151:300], ylab="Protein intensities [log2]")
boxplot(rpc[,301:ncol(rpc)], ylab="Protein intensities [log2]")

# Density plots after scaling                                                   
limma::plotDensities(rpc[,1:150], legend=FALSE)
limma::plotDensities(rpc[,151:300], legend=FALSE)
limma::plotDensities(rpc[,301:ncol(rpc)], legend=FALSE)


#########################                                                       
## Removal of outliers ##
#########################

## Per sample: Is a specific protein much more abundant within a sample?        ## 0 values removed
sum(is.na(rpc))
for (sample in colnames(rpc)){
  rpc[,sample] <- IQR_filter(rpc[,sample], 3)
}
sum(is.na(rpc))

## Per protein: Is a protein unreasonably expressed compared to other samples?  ## ncol(rpc)*nrow(rpc) = 114920
sum(is.na(rpc))                                                                 # 4993
for (protein in rownames(rpc)){
  rpc[protein,] <- IQR_filter(rpc[protein,], 3)
}
sum(is.na(rpc))                                                                 # 5139 - 4998 = 141
                                                                                # (141 / (114920-4993))*100 = 0.13%         

## Proteins identified in at least 100 samples -> statistical and biomarker considerations
qdf <- rpc
qdf[!is.na(qdf)] <- 1                                                           
hist(colSums(qdf, na.rm = T),
     xlab="Number of different proteins",
     main="Distribution of identified proteins per sample")
length(rowSums(qdf, na.rm = T)[rowSums(qdf, na.rm = T) < 100])                  # 1 proteins identified
rpc <- rpc[rowSums(qdf, na.rm = T) >= 100,]                                     # 1 proteins removed, 259 remaining


#################################################                               
## Align protein dataframe with phenotype data ##
#################################################

# Check for duplicates | No duplicates
ncol(rpc) == length(unique(colnames(rpc)))

# Replace randomization ID with patient ID + study visit
for (rID in colnames(rpc)){
  pID <- pkey$`Participant ID`[pkey$`RANDOM POSITION/ RUN NUMBER` == rID]
  tp <- pkey$`Study Visit`[pkey$`RANDOM POSITION/ RUN NUMBER` == rID]
  if (tp == "Visit 2"){
    pID <- paste(pID, "_V2", sep="")}
  else{
    pID <- paste(pID, "_V4", sep="")}
  colnames(rpc)[colnames(rpc) == rID] <- pID
}

## Align phenotype and protein datasets
# Drop data from patients without protein samples
pdf <- pdf[rownames(pdf) %in% colnames(rpc),]

# Drop protein samples without phenotype data                                   # "D011_V2" "D058_V2" "D020_V2" "D008_V2" "D019_V2" dropped
colnames(rpc)[!colnames(rpc) %in% rownames(pdf)]
rpc <- rpc[,colnames(rpc) %in% rownames(pdf)]

# Confirm similarity
table(colnames(rpc) %in% rownames(pdf))
table(rownames(pdf) %in% colnames(rpc))
ncol(rpc) == nrow(pdf)

# Reorder protein samples according to phenotype dataset
rpc <- rpc[,rownames(pdf)]
table(colnames(rpc) == rownames(pdf))

## Phenotype data
# Add technical sample information to phenotype dataset                                                     
pdf$SerumID <- gsub("_", "S", rownames(pdf))
plates <- c()
for (sID in pdf$SerumID){
  plates <- c(plates, pkey$`Plate number`[pkey$`Serum ID` == sID])
}
pdf$plate <- as.factor(plates)


###########################
## Impute missing values ##
###########################
set.seed(1)
rpci <- rpc

## Filter out proteins missing in > 20% of samples
nas <- rowSums(is.na(rpci))
rpci <- rpci[nas <= ncol(rpci)/5,]

## Transpose dataframe and impute missing values
rpci_T <- as.data.frame(t(rpci))

## Change protein names for imputation algorithm, generate key
cnames <- colnames(rpci_T)
colnames(rpci_T) <- paste("p", 1:ncol(rpci_T), sep="")
namekey <- data.frame(Original = cnames, Permuted = colnames(rpci_T))

## Setup for prediction; select 100 absolutely highest correlating predictors per protein
init=mice(rpci_T, maxit=0)
predM <- init$predictorMatrix
init=NULL

corM <- cor(rpci_T, method = "pearson", use="pairwise.complete.obs")
for (depvar in rownames(predM)){
  pred <- names(sort(abs(corM[depvar,]), decreasing=T)[2:101])
  predM[depvar, !colnames(predM) %in% pred] <- 0
}

# Impute missing values using a multithreaded version of MICE                   
rpc_Ti <- mice.par(rpci_T, predictorMatrix = predM, 
                   m = 10, maxit=50, seed = 1,
                   nnodes = 5)


#####################
## Save dataframes ##
#####################
save(rpc, rpc_Ti, pdf, namekey,
     file = "normalized_protein_intensities_p_DIA_NN.RData")

