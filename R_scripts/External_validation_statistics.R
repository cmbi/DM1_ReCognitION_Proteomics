
## Script generated on 19/04/24 by D. van As
## Purpose: Statistical analysis of of RUMC Protein Assay of external Canadian DM1 cohort
## Last update on: 11/04/25


####################
## Load libraries ##
####################
library("readxl")  #V 1.4.0


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


################
## Load files ##
################

## ELISA based ITIH3 quantification
mm <- read_excel("20231024_ITIH3_ DM1_individual values.xlsx", sheet=1)
mm <- as.data.frame(mm)

## External phenotype data (pdf)
pdf <- read_excel("Phenotype_ITIH3_DM1total.xlsx", sheet=1)
pdf <- as.data.frame(pdf)

## Log2 transformed, median centered and filtered protein intensities (All proteins = all_p_avg)
load(file = "ext_val_normalized_proteins.RDATA")


##############################################################
## Filter phenotype data for relevant measures and outliers ##
##############################################################

## Set patient ID's as rownames in phenotype dataframe
rownames(pdf) <- pdf$IDPATIENT

## Subset relevant t2 (2011) measures & simplify names
pdf <- pdf[,c(1, 4:9)]
colnames(pdf) <- c("IDPATIENT","Sex", "Age", "CTG", "SMWT",
                   "gripD", "gripG")

## Check & correct variable types
lapply(pdf, class)
pdf$Sex <- as.factor(pdf$Sex)

# Change SMWT scores from feet to meter
pdf$SMWT <- round(pdf$SMWT*0.3048,0)

## Outlier screening
# Screen all outcome measures based on 3x IQR
cdf <- pdf #Generate backup
sum(is.na(cdf))
for(OM in colnames(cdf)[3:7]){
  cdf[,OM] <- IQR_filter(cdf[,OM], 3)
}
sum(is.na(cdf))

# Manually check flagged outliers
apply(pdf, 2, function(x) length(which(!is.na(x)))) - apply(cdf, 2, function(x) length(which(!is.na(x))))

plot(sort(pdf$gripD))  ## 1 outlier removed
pdf$gripD[pdf$gripD > 50] <- NA
plot(sort(pdf$gripG))  ## 1 outlier removed
pdf$gripG[pdf$gripG > 50] <- NA


#######################################################
## Merge Canadian phenotype data with proteomic data ##
#######################################################

## Transverse protein dataframe
pr_ca_t <- as.data.frame(t(all_p_avg))

## Merge phenotype and protein dataframes
table(rownames(pr_ca_t) %in% rownames(pdf))
table(rownames(pdf) %in% rownames(pr_ca_t))                                     # 1 sample missing with ID 619
rownames(pdf)[!rownames(pdf)%in%rownames(pr_ca_t)]
pdf <- pdf[!rownames(pdf)=="619",]

pr_ca_t <- pr_ca_t[order(rownames(pr_ca_t)),]
pdf <- pdf[order(rownames(pdf)),]
table((rownames(pdf)==rownames(pr_ca_t)))
table(rownames(pr_ca_t)==rownames(pdf))

cdf <- merge(pdf, pr_ca_t, by="row.names")
rownames(cdf) <- cdf$Row.names
cdf$Row.names <- NULL


#######################################
## CTG-repeat association validation ##
#######################################

## Fit linear model for each proteingroup
fit_list <- list()
for (protein in colnames(cdf)[8:ncol(cdf)]){
  fit <- lm(formula = cdf[,protein] ~ Sex + CTG, data=cdf)
  fit_list[[protein]] <- fit
}

## Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,4)]
}
CTG_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][3,]})))

## Add FDR correction & specify column names
CTG_res_df$FDR <- p.adjust(CTG_res_df$`Pr(>|t|)`, method='fdr')
colnames(CTG_res_df) <- paste("CTG_", colnames(CTG_res_df), sep="" )
rownames(CTG_res_df) <- names(fit_coefficients)

## Add correlation coefficients
cor_vec <- c()
for (protein in rownames(CTG_res_df)){
  cor_vec <- c(cor_vec, cor(cdf[,protein], cdf[,"CTG"], use="pairwise.complete.obs", method="pearson"))
}
CTG_res_df$pcor <- cor_vec

## Quantitative analysis                                                  
sum(CTG_res_df$`CTG_Pr(>|t|)` < 0.05)                                           # 31 hits
sum(CTG_res_df$CTG_FDR < 0.05)                                                  # 0 hits
hist(CTG_res_df$`CTG_Pr(>|t|)`,
     main="P-value distribution for CTG repeat associations",
     xlab="")


#################################
## 6MWT association validation ##
#################################

# Fit linear model for each protein
fit_list <- list()
for (protein in colnames(cdf)[8:ncol(cdf)]){
  fit <- lm(formula = cdf[,protein] ~ Sex + SMWT, data=cdf)
  fit_list[[protein]] <- fit
}

# Obtain coefficient estimates and p-values
fit_coefficients <- list()
for (fit in names(fit_list)){
  fit_coefficients[[fit]] <- summary(fit_list[[fit]])[["coefficients"]][,c(1,2,4)]
}
SMWT_res_df <- as.data.frame(do.call(rbind, lapply(names(fit_coefficients), function(x){fit_coefficients[[x]][3,]})))

# Add FDR correction & specify column names
SMWT_res_df$FDR <- p.adjust(SMWT_res_df$`Pr(>|t|)`, method='fdr')
colnames(SMWT_res_df) <- paste("SMWT_", colnames(SMWT_res_df), sep="" )
rownames(SMWT_res_df) <- names(fit_coefficients)

# Add correlation coefficients
cor_vec <- c()
for (protein in rownames(SMWT_res_df)){
  cor_vec <- c(cor_vec, cor(cdf[,protein], cdf[,"SMWT"], use="pairwise.complete.obs", method="pearson"))
}
SMWT_res_df$pcor <- cor_vec

## Check for significant hits                                                   
sum(SMWT_res_df$`SMWT_Pr(>|t|)` < 0.05)                                         # 45 hits
sum(SMWT_res_df$SMWT_FDR < 0.05)                                                # 11 hits
hist(SMWT_res_df$`SMWT_Pr(>|t|)`,
     main="P-value distribution for 6MWT repeat associations",
     xlab="")


######################
## ITIH3 validation ##
######################

## Update dataframe structure, column names and change variable types
colnames(mm) <- mm[1,]
mm <- mm[-1,]
colnames(mm)[1] <- "IDPATIENT"
colnames(mm)[2] <- "ITIH3"

mm$IDPATIENT <- as.numeric(mm$IDPATIENT)
rownames(mm) <- mm$IDPATIENT
mm$age <- as.numeric(mm$age)
mm$ITIH3 <- as.numeric(mm$ITIH3)
mm$ITIH3log2 <- log2(mm$ITIH3)
mm$CTG <- round(as.numeric(mm$CTG))

## Merge datasets by PatientID
mm <- mm[!mm$IDPATIENT=="619",]
mm <- mm[order(mm$IDPATIENT),]
pdf <- pdf[order(pdf$IDPATIENT),]
table(pdf$IDPATIENT == mm$IDPATIENT)
table(mm$CTG == pdf$CTG) # correct, but rounding differences. Retain original Canadian data.
mm$CTG <- NULL
table(mm$age == pdf$Age)

c <- merge(mm, pdf, by="IDPATIENT")
rownames(c) <- c[,1]
c$Sex <- as.factor(c$Sex)

## Generate linear models based on log2ITIH3
lm_smwt_log2 <- lm(ITIH3log2 ~ Sex + SMWT, data=c)
lm_gripD_log2 <- lm(ITIH3log2 ~ Sex + gripD, data=c)
lm_gripG_log2 <- lm(ITIH3log2 ~ Sex + gripG, data=c)


##################
## Save results ##
##################

## Save results
save(c, cdf, CTG_res_df, SMWT_res_df, 
     file = "all_p_ext_val_res.RDATA")

## Generate excel output CTG association results
write.csv(CTG_res_df, file="ext_val_CTG_hits.csv")

## Generate excel output of SMWT association results
write.csv(SMWT_res_df, file="ext_val_SMWT_hits.csv")


