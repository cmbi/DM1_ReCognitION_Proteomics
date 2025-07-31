
## Script generated on 28/02/23 by D van As
## Purpose of script: ReCognitION DIA-pooled DIA-NN peptide sample pre-filter analysis
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

## Raw Peptide counts
rpc <- read.table("merged_speclib.pr_matrix.tsv", sep = "\t", header=T)

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

## Remove peptides associated with spike-in bovine-albumin and pig-trypsin PGs
any(grepl("P00764", rpc$Protein.Group))
any(grepl("TRYP_PIG", rpc$Protein.Group))
rpc[grepl("TRYP_PIG", rpc$Protein.Group),1:5]
rpc <- rpc[!grepl("TRYP_PIG", rpc$Protein.Group),]

any(grepl("P02769", rpc$Protein.Group))
any(grepl("ALBU_BOVIN", rpc$Protein.Group))
rpc[grepl("ALBU_BOVIN", rpc$Protein.Group),1:5]
rpc <- rpc[!grepl("ALBU_BOVIN", rpc$Protein.Group),]

## Simplify sample names to randomization numbers
colnames(rpc)[11:ncol(rpc)] <- gsub(".*iRT_", "", colnames(rpc)[11:ncol(rpc)])
colnames(rpc)[colnames(rpc) == "2G6_20210817142613.raw"] <- "2G6.raw"
colnames(rpc)[colnames(rpc) == "1D1_20210820165914.raw"] <- "1D1.raw"
colnames(rpc)[colnames(rpc) == "3F12_20210830074845.raw"] <- "3F12.raw"
colnames(rpc)[11:ncol(rpc)] <- gsub(".raw", "", colnames(rpc)[11:ncol(rpc)])

## Aggregation of peptides with different post-transcriptional modifications
rpc_c <- data.frame(matrix(ncol = ncol(rpc), nrow=length(unique(rpc$Stripped.Sequence))))
colnames(rpc_c) <- colnames(rpc)
rownames(rpc_c) <- unique(rpc$Stripped.Sequence)

for (peptide in rownames(rpc_c)){
  rpc_c[peptide, 11:ncol(rpc_c)] <- colSums(rpc[rpc$Stripped.Sequence == peptide, 
                                               11:ncol(rpc)], na.rm = T)
}
rpc <- rpc_c

## Drop non-relevant columns & change 0's to NA's
rpc <- rpc[,-c(1:10)]
rpc[rpc == 0] <- NA

## Set intensities of < 10.000 to NA (unreliable)
rpc[rpc < 10000] <- NA


#########################                    
## Log2 transformation ##
#########################

## Confirm that all rpcs are numeric
table(apply(rpc, 2, function(x) is.numeric(x)))

## Log2(x) transformation of intensities
rpc <- log2(rpc)


################################################################                 
## Removal of faulty samples with relatively low total counts ##
################################################################

## Total number of peptides per sample -> normal distribution
hist(colSums(rpc, na.rm=T), n = 50, 
     xlab="Total peptide abundancy per sample [Log2]",
     main="Distribution of total peptide abundancy")

# Remove samples based on sample total counts [presumed technical error]        ## 10 samples removed
lowS <- sort(colSums(rpc, na.rm=T))
lowS <- IQR_filter(lowS, 3)
rpc[,names(lowS[is.na(lowS)])] <- NULL


########################################                                        
## Median based sample normalization  ##
########################################

# Box plots before scaling
boxplot(rpc[,1:150], ylab="Peptide intensities [log2]")
boxplot(rpc[,151:300], ylab="Peptide intensities [log2]")
boxplot(rpc[,301:ncol(rpc)], ylab="Peptide intensities [log2]")

# Equal Median Normalization
rpc <- equalMedianNormalization(rpc) 

# Box plots after scaling
boxplot(rpc[,1:150], ylab="Peptide intensities [log2]")                         
boxplot(rpc[,151:300], ylab="Peptide intensities [log2]")
boxplot(rpc[,301:ncol(rpc)], ylab="Peptide intensities [log2]")

# Density plots after scaling                                                   
limma::plotDensities(rpc[,1:150], legend=FALSE)
limma::plotDensities(rpc[,151:300], legend=FALSE)
limma::plotDensities(rpc[,301:ncol(rpc)], legend=FALSE)


#########################                                                        
## Removal of outliers ##
#########################

## Per sample: Is a specific peptide much more abundant within a sample?        ## 0 values removed
sum(is.na(rpc))
for (sample in colnames(rpc)){
  rpc[,sample] <- IQR_filter(rpc[,sample], 3)
}
sum(is.na(rpc))

## Per peptide: Is a peptide unreasonably expressed compared to other samples?  # ncol(rpc) * nrow(rpc) = 1191190
sum(is.na(rpc))                                                                 # 176975
for (peptide in rownames(rpc)){
  rpc[peptide,] <- IQR_filter(rpc[peptide,], 3)
}
sum(is.na(rpc))                                                                 # 179891 - 176975 = 2916
                                                                                # (2916 / (1191190-176975)) * 100 = 0.29%  

## Number of unique peptides per sample | Peptides identified in at least 100 samples -> statistical and biomarker considerations
qdf <- rpc
qdf[!is.na(qdf)] <- 1                                                           
hist(colSums(qdf, na.rm = T),
     xlab="Number of different peptides",
     main="Distribution of identified peptides per sample")
length(rowSums(qdf, na.rm = T)[rowSums(qdf, na.rm = T) < 100])                  # 25 peptides identified
rpc <- rpc[rowSums(qdf, na.rm = T) >= 100,]                                     # 25 peptides removed, 2670 remaining


##############################################                                  
## Variance vs abundance plot after scaling ##
##############################################

pep_vars <- c()
pep_means <- c()
for (peptide in 1:nrow(rpc)){
  pep <- as.numeric(as.vector(rpc[peptide,]))
  pep_mean <- mean(pep, na.rm = T)
  pep_means <- c(pep_means, pep_mean)
  pep_var <- var(pep, na.rm = T)
  pep_vars <- c(pep_vars, pep_var)
}
plot(pep_means, pep_vars,
     main="Peptide mean-variance plot, filtered & scaled",
     xlab="Mean peptide intensity log2(x)",
     ylab="Peptide instensity variance log2(x)")


#################################################                               
## Align peptide dataframe with phenotype data ##
#################################################

## Check for duplicates | No duplicates
ncol(rpc) == length(unique(colnames(rpc)))

## Store all randomization ID's
ID_vec <- colnames(rpc)

## Replace randomization ID with patient ID + study visit
for (rID in colnames(rpc)){
  pID <- pkey$`Participant ID`[pkey$`RANDOM POSITION/ RUN NUMBER` == rID]
  tp <- pkey$`Study Visit`[pkey$`RANDOM POSITION/ RUN NUMBER` == rID]
  if (tp == "Visit 2"){
    pID <- paste(pID, "_V2", sep="")}
  else{
    pID <- paste(pID, "_V4", sep="")}
  colnames(rpc)[colnames(rpc) == rID] <- pID
}

## Align phenotype and peptide datasets
# Drop data from patients without peptide samples                               
pdf <- pdf[rownames(pdf) %in% colnames(rpc),]                                   

# Drop peptide samples without phenotype data                                   # "D011_V2" "D058_V2" "D020_V2" "D008_V2" "D019_V2"
colnames(rpc)[!colnames(rpc) %in% rownames(pdf)]
rpc <- rpc[,colnames(rpc) %in% rownames(pdf)]

# Confirm similarity
table(colnames(rpc) %in% rownames(pdf))
table(rownames(pdf) %in% colnames(rpc))
ncol(rpc) == nrow(pdf)

# Reorder peptide samples according to phenotype dataset
rpc <- rpc[,rownames(pdf)]
table(colnames(rpc) == rownames(pdf))

# Add technical sample information to phenotype dataset                                                     
pdf$SerumID <- gsub("_", "S", rownames(pdf))
plates <- c()
for (sID in pdf$SerumID){
  plates <- c(plates, pkey$`Plate number`[pkey$`Serum ID` == sID])
}
pdf$plate <- as.factor(plates)


#####################
## Save dataframes ##
#####################

save(rpc, pdf, ID_vec,
     file = "normalized_peptide_intensities_p_DIA_NN.RData")

