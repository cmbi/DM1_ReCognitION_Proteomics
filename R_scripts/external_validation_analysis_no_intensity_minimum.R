
## Script generated on 19/4/24 by D van As, last update on: 02/06/25
## Purpose: Pre-filtering of raw intensities of RUMC Protein Assay
## without intensity minimum of 10.000 - for visualization purposes
## without outlier filtering
## including german DM1 samples and controls

####################
## Load libraries ##
####################

library("arrow")         #V 7.0.0
library("DEqMS")         #V 1.8.0

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
## Load relevant dataset ##
###########################

## All quantified proteins (all_p)
all_p <- read_parquet("summary-results.proteinpgmaxlfq.parquet")
all_p <- as.data.frame(all_p)


#########################
## Data pre-processing ##
#########################

## Remove iRT-Kit spike-in P123456 
any(grepl("P123456", all_p$protein_accession))
all_p[grepl("P123456", all_p$protein_accession),1:5]
all_p <- all_p[all_p$protein_accession != "P123456",]

## Set missing values (0) to NA
all_p[all_p == 0] <- NA

## Set protein names as row names & non-relevant columns
rownames(all_p) <- all_p$protein_accession
all_p <- all_p[,-c(1:5)]

## Reduce column names to Patient IDs
colnames(all_p) <- gsub("\\_T2.*", "", colnames(all_p))
colnames(all_p) <- gsub("\\_S.*", "", colnames(all_p))
colnames(all_p) <- gsub("20240320_", "", colnames(all_p))
colnames(all_p) <- gsub("Control_", "", colnames(all_p))
colnames(all_p) <- gsub("DM_", "", colnames(all_p))
colnames(all_p) <- gsub(" ", "", colnames(all_p))

## Split technical replicates
all_p_r1 <- all_p[,c(TRUE,FALSE)]
all_p_r2 <- all_p[,c(FALSE,TRUE)]
table(colnames(all_p_r1) == colnames(all_p_r2))
corM <- cor(all_p_r1, all_p_r2, use = "pairwise.complete.obs")
hist(diag(corM), main = "Pearson rho of technical replicates")

## Calculate mean of technical replicates -> more robust than adding for missing values
all_p_avg <- matrix(nrow=nrow(all_p_r1), ncol=0)
for (patient in colnames(all_p_r1)){
  x <- cbind(all_p_r1[,patient], all_p_r2[,patient])
  x <- cbind(x, rowMeans(x, na.rm = T))
  x[,3][is.nan(x[,3])] <- NA
  all_p_avg <- cbind(all_p_avg, x[,3])
}
all_p_avg <- as.data.frame(all_p_avg)
colnames(all_p_avg) <- colnames(all_p_r1)
rownames(all_p_avg) <- rownames(all_p_r1)


#########################################                    
## Log2 transformation of protein data ##
#########################################

## Check if all rpcs are numeric | All intensities stored as numeric
table(apply(all_p_avg, 2, function(x) is.numeric(x)))

## Log2(x) transformation of intensities
all_p_avg <- log2(all_p_avg)


########################################                                        
## Median based sample normalization  ##
########################################

# Box plots before scaling
boxplot(all_p_avg, xlab="", ylab="Protein intensities [log2]")

# Equal Median Normalization
all_p_avg <- equalMedianNormalization(all_p_avg) 

# Box plots after scaling
boxplot(all_p_avg, xlab="", ylab="Protein intensities [log2]")


###########################
## Save complete dataset ##
###########################
save(all_p_avg,
     file = "ext_val_normalized_proteins_no_minimum.RDATA")

