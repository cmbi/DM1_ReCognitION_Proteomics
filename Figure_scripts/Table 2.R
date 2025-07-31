
## Script generated on 30/7/24 by D van As, last update on: 02/06/25
## Table 2: Combined validated CTG and 6MWT associations

################
## Load files ##
################

## ReCognitION statistical modeling results
load(file = "DIA-NN_p_Protein_statistics.RDATA")
rec_CTG <- Mode_res_df
rec_SMWT <- outcome_coef[["SMWT"]]
rec_DM1ActivC <- outcome_coef[["DM1ActivC"]]
rec_MeanENMO <- outcome_coef[["MeanENMO"]]
rec_M5ENMO <- outcome_coef[["M5ENMO"]]
rm(Vis_Treat_res_df, Vis_Get_res_df, Mode_res_df, outcome_coef, coef_table)

## External validation statistical modeling results
load(file = "all_p_ext_val_res.RDATA")
ext_CTG <- CTG_res_df
ext_SMWT <- SMWT_res_df
rm(CTG_res_df, SMWT_res_df, c, cdf)


###############################################################
## Table 2: Combined validated CTG and 6MWT associations ##
###############################################################

#### CTG-repeat associations
## Identify overlapping proteins 

# Obtain proteins from significant (FDR based) ReCognitION protein groups
rec_sig_proteins <- c()
rec_sig_PG <- rownames(rec_CTG[rec_CTG$Mode_FDR < 0.05,])
for (PG in rec_sig_PG){
  rec_sig_proteins <- c(rec_sig_proteins, strsplit(PG, ";")[[1]])
}
rec_sig_proteins <- unique(rec_sig_proteins)

# Obtain proteins from external validation nominal significant protein groups
ext_sig_proteins <- c()
ext_sig_PG <- rownames(ext_CTG[ext_CTG$`CTG_Pr(>|t|)`<0.05,])
for (PG in ext_sig_PG){
  ext_sig_proteins <- c(ext_sig_proteins, strsplit(PG, ";")[[1]])
}
ext_sig_proteins <- unique(ext_sig_proteins)

# Intersect
CTG_hits_int <- intersect(rec_sig_proteins, ext_sig_proteins)

## For each protein, identify the most significant protein-group in OPTIMISTIC
CTG_val <- data.frame(matrix(ncol = ncol(rec_CTG)+1, nrow=0))
colnames(CTG_val) <- c(colnames(rec_CTG), "Protein_group")

for (protein in CTG_hits_int){
  df <- rec_CTG[grepl(protein, rownames(rec_CTG)),]
  df <- df[which.min(df$`Mode_Pr(>|t|)`),]
  df$Protein_group <- rownames(df)
  CTG_val <- rbind(CTG_val, df)
}

# Reduce identified protein groups to unique hits
CTG_val <- unique(CTG_val)

## For each selected OPTIMISTIC protein group, find the best match in the external validation set

CTG_val$ext_protein_group <- NA
CTG_val$ext_p.value <- NA
CTG_val$ext_pcorr <- NA

for (rPG in CTG_val$Protein_group){
  rec_prot <- strsplit(rPG, ";")[[1]]
  
  max_overlap <- 0  # Keeps track of maximum protein overlap
  max_prot <- 0     # Keeps track of the number of proteins in the external protein(group)
  min_sig <- 1      # Keeps track of p-values
  match <- c()      # Stores best match(es)
  
  ## Obtain all (partially) matching protein groups in significant external data
  for (ePG in ext_sig_PG){
    sig <- ext_CTG$`CTG_Pr(>|t|)`[rownames(ext_CTG)==ePG]
    ext_prot <- strsplit(ePG, ";")[[1]]
    n_prot <- length(ext_prot)
    overlap <- length(intersect(rec_prot, ext_prot))
    
    # Replace match if higher overlap is found
    if (overlap > max_overlap){
      match <- ePG
      max_overlap <- overlap
      max_prot <- n_prot
      min_sig <- sig
      
    # Replace match if overlap is equal but external protein group is smaller
    } else if (overlap > 0 & overlap == max_overlap & n_prot < max_prot){
      match <- ePG
      max_prot <- n_prot
      min_sig <- sig
      
    # Replace match if overlap and pg size is equal, but more significant  
    } else if (overlap > 0 & overlap == max_overlap & n_prot == max_prot & sig < min_sig){
      match <- ePG
      min_sig <- sig
      
    # Add match if new result is equal to previous best result
    } else if (overlap > 0 & overlap == max_overlap & n_prot == max_prot & sig == min_sig){
      match <- c(match, ePG)
    }
  }
  
  ## Select optimal match & Add to results
  if (length(match) == 0){
    CTG_val[CTG_val$Protein_group==rPG, "ext_protein_group"] <- "No match"
  }
  else if (length(match) == 1){
    CTG_val[CTG_val$Protein_group==rPG, "ext_protein_group"] <- match
    CTG_val[CTG_val$Protein_group==rPG, "ext_p.value"] <- ext_CTG$`CTG_Pr(>|t|)`[rownames(ext_CTG)==match]
    CTG_val[CTG_val$Protein_group==rPG, "ext_pcorr"] <- ext_CTG$pcor[rownames(ext_CTG)==match]
  }
  else if (length(match) > 1) {
    CTG_val[CTG_val$Protein_group==rPG, "ext_protein_group"] <- "Multiple matches"
  } 
}  

## Drop & reorder columns
CTG_val$Mode_Estimate <- NULL
CTG_val$`Mode_Std. Error` <- NULL
CTG_val$`Mode_Pr(>|t|)` <- NULL
CTG_val <- CTG_val[,c(3, 1, 2, 4, 5, 6)]
CTG_val <- CTG_val[order(abs(CTG_val$pcor), decreasing = T),]

## Round values to 4 digits
CTG_val[,c(2, 3, 5, 6)] <-
  apply(CTG_val[,c(2, 3, 5, 6)], 2, function(x){
    round(x, digits = 4)
  })

## Change column names
rownames(CTG_val) <- NULL
colnames(CTG_val) <- c("Rec_PG", "Rec_FDR", "Rec_pcor",
                       "Ext_PG", "Ext_pval", "Ext_pcor")

#### 6MWT Associations

# Obtain proteins from significant (FDR based) ReCognitION protein groups
rec_sig_proteins <- c()
rec_SMWT <- rec_SMWT[rec_SMWT$SMWT_FDR < 0.05, ]
for (PG in rec_sig_PG){
  rec_sig_proteins <- c(rec_sig_proteins, strsplit(PG, ";")[[1]])
}
rec_sig_proteins <- unique(rec_sig_proteins)

# Obtain proteins from external validation nominal significant protein groups
ext_sig_proteins <- c()
ext_sig_PG <- rownames(ext_SMWT[ext_SMWT$`SMWT_Pr(>|t|)` < 0.05,])
for (PG in ext_sig_PG){
  ext_sig_proteins <- c(ext_sig_proteins, strsplit(PG, ";")[[1]])
}
ext_sig_proteins <- unique(ext_sig_proteins)

# Intersect
SMWT_hits_int <- intersect(rec_sig_proteins, ext_sig_proteins)                  # N=19

## For each protein, identify the most significant protein-group in OPTIMISTIC
SMWT_val <- data.frame(matrix(ncol = ncol(rec_SMWT)+1, nrow=0))
colnames(SMWT_val) <- c(colnames(rec_SMWT), "Protein_group")

for (protein in SMWT_hits_int){
  df <- rec_SMWT[grepl(protein, rownames(rec_SMWT)),]
  df <- df[which.min(df$`SMWT_Pr(>|t|)`),]
  df$Protein_group <- rownames(df)
  SMWT_val <- rbind(SMWT_val, df)
}

# Reduce identified protein groups to unique hits
SMWT_val <- unique(SMWT_val)

## For each selected OPTIMISTIC protein group, find the best match in the external validation set
SMWT_val$ext_protein_group <- NA
SMWT_val$ext_p.value <- NA
SMWT_val$ext_pcorr <- NA

for (rPG in SMWT_val$Protein_group){
  rec_prot <- strsplit(rPG, ";")[[1]]
  
  max_overlap <- 0  # Keeps track of maximum protein overlap
  max_prot <- 0     # Keeps track of the number of proteins in the external protein(group)
  min_sig <- 1      # Keeps track of p-values
  match <- c()      # Stores best match(es)
  
  ## Obtain all (partially) matching protein groups in significant external data
  for (ePG in ext_sig_PG){
    sig <- ext_SMWT$`SMWT_Pr(>|t|)`[rownames(ext_SMWT)==ePG]
    ext_prot <- strsplit(ePG, ";")[[1]]
    n_prot <- length(ext_prot)
    overlap <- length(intersect(rec_prot, ext_prot))
    
    # Replace match if higher overlap is found
    if (overlap > max_overlap){
      match <- ePG
      max_overlap <- overlap
      max_prot <- n_prot
      min_sig <- sig
      
      # Replace match if overlap is equal but external protein group is smaller
    } else if (overlap > 0 & overlap == max_overlap & n_prot < max_prot){
      match <- ePG
      max_prot <- n_prot
      min_sig <- sig
      
      # Replace match if overlap and pg size is equal, but more significant  
    } else if (overlap > 0 & overlap == max_overlap & n_prot == max_prot & sig < min_sig){
      match <- ePG
      min_sig <- sig
      
      # Add match if new result is equal to previous best result
    } else if (overlap > 0 & overlap == max_overlap & n_prot == max_prot & sig == min_sig){
      match <- c(match, ePG)
    }
  }
  
  ## Select optimal match & Add to results
  if (length(match) == 0){
    SMWT_val[SMWT_val$Protein_group==rPG, "ext_protein_group"] <- "No match"
  }
  else if (length(match) == 1){
    SMWT_val[SMWT_val$Protein_group==rPG, "ext_protein_group"] <- match
    SMWT_val[SMWT_val$Protein_group==rPG, "ext_p.value"] <- ext_SMWT$`SMWT_Pr(>|t|)`[rownames(ext_SMWT)==match]
    SMWT_val[SMWT_val$Protein_group==rPG, "ext_pcorr"] <- ext_SMWT$pcor[rownames(ext_SMWT)==match]
  }
  else if (length(match) > 1) {
    SMWT_val[SMWT_val$Protein_group==rPG, "ext_protein_group"] <- "Multiple matches"
  } 
}  

## Drop & reorder columns
SMWT_val$SMWT_Estimate <- NULL
SMWT_val$`SMWT_Std. Error` <- NULL
SMWT_val$`SMWT_Pr(>|t|)` <- NULL
SMWT_val <- SMWT_val[,c(3, 1, 2, 4, 5, 6)]
SMWT_val <- SMWT_val[order(abs(SMWT_val$pcor), decreasing = T),]

## Round values to 4 digits
SMWT_val[,c(2, 3, 5, 6)] <-
  apply(SMWT_val[,c(2, 3, 5, 6)], 2, function(x){
    round(x, digits = 4)
  })

## Change column names
rownames(SMWT_val) <- NULL
colnames(SMWT_val) <- c("Rec_PG", "Rec_FDR", "Rec_pcor",
                       "Ext_PG", "Ext_pval", "Ext_pcor")

## Combine results into one table T2
CTG_val$OM <- "CTG"
SMWT_val$OM <- "6MWT"
table(colnames(CTG_val) == colnames(SMWT_val))
T2 <- rbind(CTG_val, SMWT_val)
T2 <- T2[,c(7, 1:6)]

## Exclude hits with opposing correlation coefficients [-> not validated]
RN <- c()
for (i in 1:nrow(T2)){
  Rec_coef <- T2[i, "Rec_pcor"]
  Ext_coef <- T2[i, "Ext_pcor"]
  if (Rec_coef > 0 & Ext_coef < 0 | Rec_coef < 0 & Ext_coef > 0)
    RN <- c(RN, i)
}
T2 <- T2[-RN,]
rownames(T2) <- NULL

## Save results
write.csv(T2, file="Table2.csv")
save(T2, 
     file = "Table2.RDATA")


##################################################################
## Study overlap in (capacity for activity) associated findings ##
##################################################################

## Subset significant hits
rec_DM1ActivC <- rec_DM1ActivC[rec_DM1ActivC$DM1ActivC_FDR < 0.05, ]
rec_MeanENMO <- rec_MeanENMO[rec_MeanENMO$MeanENMO_FDR < 0.05, ]
rec_M5ENMO <- rec_M5ENMO[rec_M5ENMO$M5ENMO_FDR < 0.05,]

## Intersect results
nrow(rec_DM1ActivC)
length(intersect(rownames(rec_SMWT), rownames(rec_DM1ActivC)))
intersect(rownames(rec_DM1ActivC), T2$Rec_PG[T2$OM=="6MWT"])
setdiff(T2$Rec_PG[T2$OM=="6MWT"], rownames(rec_DM1ActivC))

nrow(rec_MeanENMO)
length(intersect(rownames(rec_SMWT), rownames(rec_MeanENMO)))
intersect(rownames(rec_MeanENMO), T2$Rec_PG[T2$OM=="6MWT"])
setdiff(T2$Rec_PG[T2$OM=="6MWT"], rownames(rec_MeanENMO))

nrow(rec_M5ENMO)
length(intersect(rownames(rec_SMWT), rownames(rec_M5ENMO)))
intersect(rownames(rec_M5ENMO), T2$Rec_PG[T2$OM=="6MWT"])
setdiff(T2$Rec_PG[T2$OM=="6MWT"], rownames(rec_M5ENMO))

