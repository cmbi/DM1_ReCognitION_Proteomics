
## Script generated on 27/7/24 by D van As, last update on: 04/06/25
## Figure 1: OPTIMISTIC protein dataset PCA
## SFigure 3: Proteomics technical analyses
## SFigure 4: OPTIMISTIC peptide dataset PCA

####################
## Load libraries ##
####################
library(ggplot2)       #V 3.3.6
library(cowplot)       #V 1.1.1
library(ggcorrplot)    #V 0.1.4
library(lme4)          #V 1.1.29
library(MuMIn)         #V 1.46.0
library(RColorBrewer)  #V 1.1.3


################
## Load files ##
################

## ReCognitION peptide data
load(file = "normalized_peptide_intensities_p_DIA_NN.RData")
pep_pdf <- pdf
pep_rpc <- rpc
rm(pdf, rpc)

## ReCognitION proteomics data
load(file = "normalized_protein_intensities_p_DIA_NN.RData")

## External validation results
load(file = "ext_val_normalized_proteins.RDATA")

## Set ggplot2 theme
ptheme <-   theme(
  aspect.ratio = 1,
  legend.title = element_text(size=14, face="bold"),
  legend.text = element_text(face="bold",color = "black", size = 12),
  legend.key = element_rect(fill="white"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  panel.background = element_rect(fill="white"),
  panel.grid.major.x = element_line(size = 0.25, color = "grey"),
  panel.grid.major.y = element_line(size = 0.25, color = "grey"),
  axis.text.x = element_text(color = "black", size = 16),
  axis.text.y = element_text(color = "black", size = 16),
  axis.title.y = element_text(face="bold", size = 16, color = "black"),
  axis.title.x = element_text(face="bold", size = 16, color = "black"),
  plot.margin = margin(t=5,r=5,b=5,l=5,unit="pt")
)

## Set color scheme
col <- colorRampPalette(c("yellow", "red"))
col <- col(20)


##################################
## Fig 1 PCA on protein dataset ##
##################################

## PCA Analysis
# Exclude patients with substantial missing proteins [> 25%] & select complete cases
qdf <- rpc
qdf[!is.na(qdf)] <- 1 
quantile(colSums(qdf, na.rm = T))

rpc_c <- rpc[,colSums(qdf, na.rm = T) > 245]                                    # 315 patients remaining
rpc_c <- rpc_c[complete.cases(rpc_c),]                                          # 173 proteins remaining

# Run on complete file
rpc_c_t <- as.data.frame(t(rpc_c)) 
pca_rpc <- prcomp(rpc_c_t, scale=T,center=T)
pca_df <- as.data.frame(pca_rpc$x)
var_explained <- data.frame(PC = paste0("PC", 1:173),
                            var_explained=(pca_rpc$sdev)^2/sum((pca_rpc$sdev)^2))

# Add phenotype data to PCA df                                                 
table(rownames(pdf) == rownames(pca_df))
pdf_pca <- pdf[rownames(pdf)%in%rownames(pca_df),]
table(rownames(pca_df) == rownames(pdf_pca))
pca_df <- merge(pca_df[,1:10], pdf_pca, by ="row.names")
rownames(pca_df) <- pca_df$Row.names
pca_df$Row.names <- NULL
pca_df$PatientID <- gsub("_V.", "", rownames(pca_df))

# Calculate R-sq values
pca_tests <- data.frame() 
for (OM in colnames(pca_df)[c(11:43, 45)]){
  for (PCA in colnames(pca_df)[1:10]){
    formula <- paste(PCA, OM, sep=" ~ ")
    formula <- paste0(formula, " + (1|PatientID)")
    fit <- suppressMessages(lmer(formula = formula, data=pca_df))
    pca_tests[PCA,OM] <- r.squaredGLMM(fit)[1]
  }
}

## Generate correlogram of marginal R-sq values
# Change variable names for plotting & order descending
colnames(pca_tests)[1] <- "AES-C"
colnames(pca_tests)[2] <- "AES-I"
colnames(pca_tests)[6] <- "CISActivity"
colnames(pca_tests)[11] <- "GET"
colnames(pca_tests)[33] <- "Study timepoint"
colnames(pca_tests)[31] <- "Study group"
colnames(pca_tests)[26] <- "SSLD"
colnames(pca_tests)[27] <- "SSLI"
colnames(pca_tests)[28] <- "SSLN"
colnames(pca_tests)[29] <- "SCWTi"
colnames(pca_tests)[25] <- "6MWT"
colnames(pca_tests)[24] <- "Study centre"
colnames(pca_tests)[23] <- "Sex"
colnames(pca_tests)[21] <- "CTG Repeat"
colnames(pca_tests)[34] <- "Well plate"
pca_tests <- pca_tests[,order(colnames(pca_tests), decreasing = T)]

A <- ggcorrplot(pca_tests)+
  scale_fill_gradient2(
    low="red",
    mid="white",
    high ="darkblue", 
    breaks=seq(0,0.3,0.1), limit=c(0,0.3))+
  labs(fill="Marginal R-squared")+
  theme(
    aspect.ratio = 2.5,
    legend.position="top",
    legend.justification = "right",
    legend.title = element_text(size=14, face="bold"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(face="bold", color = "black", size = 15),
    axis.text.y = element_text(face="bold", color = "black", size = 15),
    plot.margin = margin(t=5,r=5,b=5,l=5,unit="pt")
  )

## Scree plot
var_explained$PC <- 1:nrow(var_explained)
B <- ggplot(var_explained[1:10,], 
       aes(x=factor(PC, level=c("1", "2", "3", "4", "5",
                                "6", "7", "8", "9", "10")),
           y=var_explained, group=1))+
  geom_point(size=2)+
  geom_line()+
  xlab("Principal component") +
  ylab("Variance explained") +
  ptheme

## PC1 and 2: CTG repeat
C <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Mode)) +
  geom_point() +
  scale_color_gradientn('CTG Repeat',colors=col)+
  ptheme

## PC2 and 3: 6MWT
D <- ggplot(pca_df, aes(x=PC2, y=PC3, color=SMWT)) +
  geom_point() +
  scale_color_gradientn('6MWT', colors=col)+
  ptheme

## PC 4 and 5: Biological Sex
colnames(pca_df)[33] <- "Sex"                                                   
pca_df$Sex <- ifelse(pca_df$Sex==0, "Female", "Male")
E <- ggplot(pca_df, aes(x=PC4, y=PC5, color=Sex)) +
  geom_point() +
  scale_y_continuous(limits = c(-6,10)) +
  ptheme +
  theme(legend.position = c(0.60, 0.92),
        legend.box.just = "right",
        legend.direction = "horizontal")

## Cast subplots into one and save
p1 <- plot_grid(A,
                labels=c("A"),
                label_size = 20)

p2 <- plot_grid(B, E,
          ncol=1,
          align="v",
          axis="lr",
          labels=c("B", "D"),
          label_size = 20)

p3 <- plot_grid(C, D,
                ncol=1,
                align="v",
                axis="r",
                labels=c("C", "E"),
                label_size = 20)


Fig1 <- plot_grid(p1, p2, p3, nrow=1,
                  rel_widths = c(1, 1, 1.1))

ggsave(Fig1, 
       file ="Fig1_proteomics_pca.png", 
       units="cm", height = 20, width = 40, 
       dpi = 600,
       device = "png")


###############################
## SFig 3 Technical analyses ##
###############################

### Within vs between sample variance
## Obtain paired samples
id_temp <- gsub("_V.*", "", colnames(rpc))
paired_IDs <- table(id_temp)
paired_IDs <- paired_IDs[paired_IDs == 2]
paired_IDs <- names(paired_IDs)                                                 # 189 samples
paired_IDs <- paired_IDs[order(paired_IDs)]

paired_V2 <- paste0(paired_IDs, "_V2")
paired_V4 <- paste0(paired_IDs, "_V4")
paired_columns <- c(paired_V2, paired_V4)

## Subset paired protein samples and confirm correct order
paired_rpc <- rpc[,colnames(rpc) %in% paired_columns]
table(colnames(paired_rpc) == paired_columns)

## Generate correlation matrix
paired_rpc_cor <- cor(paired_rpc, method="pearson", use="pairwise.complete.obs") 
paired_rpc_cor <- paired_rpc_cor[1:189,190:378]

## Identify samples with highest within sample correlation                      # 143 samples
mc_ids <- c()
for (i in 1:nrow(paired_rpc_cor)){
  d <- paired_rpc_cor[i,i]
  m <- max(paired_rpc_cor[i,-i], na.rm = T)
  if (d > m){
    mc_ids <- c(mc_ids, rownames(paired_rpc_cor)[i])
  }
}

## Obtain within and between sample correlations
sc <- diag(paired_rpc_cor)

# Set self correlations to NA
for (i in 1:nrow(paired_rpc_cor)){
  paired_rpc_cor[i,i] <- NA
}
# Obtain average correlations across samples & visualize
avg_cor <- rowMeans(paired_rpc_cor, na.rm = T)

## Generate dataframe for plotting
Rho <- c(sc, avg_cor)
Label <- c(rep("Paired Rho", length(sc)),
           rep("Average Rho", length(avg_cor)))
cor_df <- as.data.frame(cbind(Rho, Label))
cor_df$Rho <- as.numeric(cor_df$Rho)

## Generate plot
A <- ggplot(data = cor_df, aes(x=Rho, fill=Label)) +
  geom_histogram(bins=50, color="black", alpha=0.6, position = "identity") +
  scale_x_continuous(breaks=seq(0.86,1,0.04), 
                     limits = c(0.85,1)) +
  xlab("Pearson Rho") +
  ylab("Frequency")+
  ptheme +
  theme(
    legend.position = c(0.30, 0.85),
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white")
  )
  

### Sample distributions over wells plates
pdf$plate <- as.numeric(as.character(pdf$plate))
pdf$Visit <- ifelse(pdf$Visit=="V2", "Baseline", "10 Months")

B <- ggplot(data = pdf, aes(x=plate, fill=Visit)) +
  geom_histogram(bins=5, color="black", alpha=0.6, position = "identity") +
  xlab("Wells plate") +
  ylab("Number of samples")+
  scale_y_continuous(breaks=seq(0,100,25), 
                     limits = c(0,115)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.92),
    legend.box.just = "right",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.key = element_rect(fill="white")
  )

### Mean variance plots
## OPTIMISTIC data
# Calculate mean and variance for each protein across all samples
prot_vars <- c()
prot_means <- c()
for (protein in 1:nrow(rpc)){
  prot <- as.numeric(as.vector(rpc[protein,]))
  prot_mean <- mean(prot, na.rm = T)
  prot_means <- c(prot_means, prot_mean)
  prot_var <- var(prot, na.rm = T)
  prot_vars <- c(prot_vars, prot_var)
}
# Generate mean/variance plot
rec_var_df <- as.data.frame(cbind(prot_means, prot_vars))

C <- ggplot(rec_var_df, aes(x=prot_means, y=prot_vars))+
  geom_point(size=1) +
  xlab("Mean PG intensity (log2)")+
  ylab("PG intensity variance (log2)")+
  ptheme 


## External validation data
# Calculate mean and variance for each protein across all samples
prot_vars <- c()
prot_means <- c()
for (protein in 1:nrow(all_p_avg)){
  prot <- as.numeric(as.vector(all_p_avg[protein,]))
  prot_mean <- mean(prot, na.rm = T)
  prot_means <- c(prot_means, prot_mean)
  prot_var <- var(prot, na.rm = T)
  prot_vars <- c(prot_vars, prot_var)
}
# Generate mean/variance plot
ext_var_df <- as.data.frame(cbind(prot_means, prot_vars))

D <- ggplot(ext_var_df, aes(x=prot_means, y=prot_vars))+
  geom_point(size=1) +
  xlab("Mean PG intensity (log2)")+
  ylab("PG intensity variance (log2)")+
  ptheme 

## Merge into one plot and save

# Arrange CTG plots in grid & save
SFig3 <- plot_grid(A, B, C, D,
          ncol=2,
          labels="AUTO",
          label_size = 20,
          align = "hv")
          
ggsave(SFig3, 
       file ="SFig3_Technical_analyses.png", 
       units="cm", height = 20, width = 24, 
       dpi = 600,
       device = "png")


###################################
## SFig 4 PCA on peptide dataset ##
###################################

## PCA analysis
# Exclude patients with substantial missing peptides [> 25%] & select complete cases
qdf <- pep_rpc
qdf[!is.na(qdf)] <- 1 
quantile(colSums(qdf, na.rm = T))

rpc_c <- pep_rpc[,colSums(qdf, na.rm = T) > 2207]                               # 327 patients remaining
rpc_c <- rpc_c[complete.cases(rpc_c),]                                          # 901 peptides remaining

# Run on complete file
rpc_c_t <- as.data.frame(t(rpc_c)) 
pca_rpc <- prcomp(rpc_c_t, scale=T,center=T)
pca_df <- as.data.frame(pca_rpc$x)
var_explained <- data.frame(PC = paste0("PC", 1:327),
                            var_explained=(pca_rpc$sdev)^2/sum((pca_rpc$sdev)^2))

# Add phenotype data to PCA df                                                 
table(rownames(pep_pdf) == rownames(pca_df))
pdf_pca <- pep_pdf[rownames(pep_pdf) %in% rownames(pca_df),]
table(rownames(pdf_pca) == rownames(pca_df))
pca_df <- merge(pca_df[,1:10], pdf_pca, by ="row.names")
rownames(pca_df) <- pca_df$Row.names
pca_df$Row.names <- NULL
pca_df$PatientID <- gsub("_V.", "", rownames(pca_df))

# Calculate R-sq values
pca_tests <- data.frame() 
for (OM in colnames(pca_df)[c(11:43, 45)]){
  for (PCA in colnames(pca_df)[1:10]){
    formula <- paste(PCA, OM, sep=" ~ ")
    formula <- paste0(formula, " + (1|PatientID)")
    fit <- suppressMessages(lmer(formula = formula, data=pca_df))
    pca_tests[PCA,OM] <- r.squaredGLMM(fit)[1]
  }
}

## Generate correlogram of marginal R-sq values
colnames(pca_tests)[1] <- "AES-C"
colnames(pca_tests)[2] <- "AES-I"
colnames(pca_tests)[6] <- "CISActivity"
colnames(pca_tests)[11] <- "GET"
colnames(pca_tests)[33] <- "Study timepoint"
colnames(pca_tests)[31] <- "Study group"
colnames(pca_tests)[26] <- "SSLD"
colnames(pca_tests)[27] <- "SSLI"
colnames(pca_tests)[28] <- "SSLN"
colnames(pca_tests)[29] <- "SCWTi"
colnames(pca_tests)[25] <- "6MWT"
colnames(pca_tests)[24] <- "Study centre"
colnames(pca_tests)[23] <- "Sex"
colnames(pca_tests)[21] <- "CTG Repeat"
colnames(pca_tests)[34] <- "Well plate"
pca_tests <- pca_tests[,order(colnames(pca_tests), decreasing = T)]

A <- ggcorrplot(pca_tests)+
  scale_fill_gradient2(
    low="red",
    mid="white",
    high ="darkblue", 
    breaks=seq(0,0.3,0.1), limit=c(0,0.3))+
  labs(fill="Marginal R-squared")+
  theme(
    aspect.ratio = 2.5,
    legend.position="top",
    legend.justification = "right",
    legend.title = element_text(size=14, face="bold"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(face="bold", color = "black", size = 15),
    axis.text.y = element_text(face="bold", color = "black", size = 15),
    plot.margin = margin(t=5,r=5,b=5,l=5,unit="pt")
  )

## Scree plot
var_explained$PC <- 1:nrow(var_explained)
B <- ggplot(var_explained[1:10,], 
            aes(x=factor(PC, level=c("1", "2", "3", "4", "5",
                                     "6", "7", "8", "9", "10")),
                y=var_explained, group=1))+
  geom_point(size=2)+
  geom_line()+
  xlab("Principal component") +
  ylab("Variance explained") +
  ptheme

## PC1 and 2: Well plate
C <- ggplot(pca_df, aes(x=PC1, y=PC2, color=plate)) +
  geom_point() +
  ptheme+
  labs(color="Well plate")

## PC4 and 5: CTG repeat
D <- ggplot(pca_df, aes(x=PC4, y=PC5, color=Mode)) +
  geom_point() +
  scale_y_continuous(limits = c(-22,22)) +
  scale_color_gradientn('CTG Repeat', colors=col)+
  ptheme

## PC8 and 9: 
colnames(pca_df)[33] <- "Sex"                                                   
pca_df$Sex <- ifelse(pca_df$Sex==0, "Female", "Male")
E <- ggplot(pca_df, aes(x=PC8, y=PC9, color=Sex)) +
  geom_point() +
  scale_y_continuous(limits = c(-12,23)) +
  ptheme +
  theme(legend.position = c(0.60, 0.92),
        legend.box.just = "right",
        legend.direction = "horizontal")


## Cast subplots into one and save
p1 <- plot_grid(A,
                labels=c("A"),
                label_size = 20)

p2 <- plot_grid(B, E,
                ncol=1,
                align="v",
                axis="lr",
                labels=c("B", "D"),
                label_size = 20)

p3 <- plot_grid(D, C,
                ncol=1,
                align="v",
                axis="r",
                labels=c("C", "E"),
                label_size = 20)

SFig4 <- plot_grid(p1, p2, p3, nrow=1,
                  rel_widths = c(1, 1, 1.1))

ggsave(SFig4, 
       file ="SFig4_pepeptides_pca.png", 
       units="cm", height = 20, width = 40, 
       dpi = 600,
       device = "png")

