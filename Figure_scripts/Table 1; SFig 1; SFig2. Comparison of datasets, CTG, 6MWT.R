
## Script generated on 27/7/24 by D van As, last update on: 09/11/24
## Table 1: Comparative features of OPTIMISTIC and Validation proteomics datasets
## SFigure 1: Differences in CTG repeat distributions between the two cohorts
## SFigure 2: Differences in 6MWT score distributions between the two cohorts


####################
## Load libraries ##
####################
library(ggplot2)  #V 3.3.6
library(MASS)     #V 7.3.53
library(cowplot)  #V 1.1.1


###############
## Functions ##
###############

## Function to apply boxcox transformation with optimal lambda
boxcox_transformation <- function(x){
  model <- boxcox(x ~ 1)
  lambda <- model$x[which.max(model$y)]
  if (lambda == 0){
    new_x <- log(x)
  }
  else{
    new_x <- (((x^lambda)-1) / lambda)
  }
  return(new_x)
}


################
## Load files ##
################

## ReCognitION results
load(file = "Generated_datasets/DIA-NN_p_Protein_statistics.RDATA")
rec_cdf <- cdf
V2_rec_cdf <- rec_cdf[rec_cdf$Visit=="V2",]
V4_rec_cdf <- rec_cdf[rec_cdf$Visit=="V4",]
rm(cdf, coef_table, Mode_res_df, outcome_coef, Vis_Get_res_df, Vis_Treat_res_df)

## External validation results
load(file = "all_p_ext_val_res.RDATA")
ext_cdf <- cdf
rm(cdf, c, CTG_res_df, SMWT_res_df)

## Set ggplot2 theme
ptheme <- theme(
  aspect.ratio = 1,
  legend.title = element_text(size=14, face="bold"),
  legend.text = element_text(color = "black", size = 12),
  legend.key = element_rect(fill="white"),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  panel.background = element_rect(fill="white"),
  panel.grid.major.x = element_line(size = 0.25, color = "grey"),
  panel.grid.major.y = element_line(size = 0.25, color = "grey"),
  plot.title = element_text(face="bold", color="black", size=16),
  axis.text.x = element_text(color = "black", size = 16),
  axis.text.y = element_text(color = "black", size = 16),
  axis.title.y = element_text(face="bold", size = 16, color = "black"),
  axis.title.x = element_text(face="bold", size = 16, color = "black"),
)

#################################
## Table 1: Cohort comparisons ##
#################################

## Obtain comparative information for Table 1
# OPTIMISTIC baseline
nrow(V2_rec_cdf)
table(V2_rec_cdf$SexCode)
min(V2_rec_cdf$Age, na.rm = T)
max(V2_rec_cdf$Age, na.rm = T)
median(V2_rec_cdf$Age, na.rm = T)
min(V2_rec_cdf$Mode, na.rm = T)
max(V2_rec_cdf$Mode, na.rm = T)
median(V2_rec_cdf$Mode, na.rm = T)
min(V2_rec_cdf$SMWT, na.rm = T)
max(V2_rec_cdf$SMWT, na.rm = T)
median(V2_rec_cdf$SMWT, na.rm = T)

# OPTIMISTIC 10 months
nrow(V4_rec_cdf)
table(V4_rec_cdf$SexCode)
min(V4_rec_cdf$Age, na.rm = T)
max(V4_rec_cdf$Age, na.rm = T)
median(V4_rec_cdf$Age, na.rm = T)
min(V4_rec_cdf$Mode, na.rm = T)
max(V4_rec_cdf$Mode, na.rm = T)
median(V4_rec_cdf$Mode, na.rm = T)
min(V4_rec_cdf$SMWT, na.rm = T)
max(V4_rec_cdf$SMWT, na.rm = T)
median(V4_rec_cdf$SMWT, na.rm = T)

# External dataset
nrow(ext_cdf)
table(ext_cdf$Sex)
min(ext_cdf$Age, na.rm = T)
max(ext_cdf$Age, na.rm = T)
median(ext_cdf$Age, na.rm = T)
length(colnames(ext_cdf)[8:ncol(ext_cdf)])
min(ext_cdf$CTG, na.rm = T)
max(ext_cdf$CTG, na.rm = T)
median(ext_cdf$CTG, na.rm = T)
min(ext_cdf$SMWT, na.rm = T)
max(ext_cdf$SMWT, na.rm = T)
median(ext_cdf$SMWT, na.rm = T)

## DM1 German cohort                                                            ## Missing information!
# See Complete Sample List 
G_DM1_ages <- c(11, 18, 33, 12, 10, 38, 20, 5, 36, 39, 41, 12, 16)
length(G_DM1_ages)
min(G_DM1_ages, na.rm = T)
max(G_DM1_ages, na.rm = T)
median(G_DM1_ages, na.rm = T)

## German control cohort
G_C_ages <- c(1, 4, 12, 4, 53, 15, 12, 28, 34, NA)
length(G_C_ages)
min(G_C_ages, na.rm = T)
max(G_C_ages, na.rm = T)
median(G_C_ages, na.rm = T)


###################################################################################
## SFig1 and 2: Comparison of (scaled) CTG and 6MWT distributions across cohorts ##
###################################################################################

## Generate scaled results
V2_rec_cdf$Mode_s <- scale(boxcox_transformation(V2_rec_cdf$Mode))[1:nrow(V2_rec_cdf)]
V2_rec_cdf$SMWT_s <- scale(boxcox_transformation(V2_rec_cdf$SMWT))[1:nrow(V2_rec_cdf)]

V4_rec_cdf$Mode_s <- scale(boxcox_transformation(V4_rec_cdf$Mode))[1:nrow(V4_rec_cdf)]
V4_rec_cdf$SMWT_s <- scale(boxcox_transformation(V4_rec_cdf$SMWT))[1:nrow(V4_rec_cdf)]

ext_cdf$CTG_s <- scale(boxcox_transformation(ext_cdf$CTG))[1:nrow(ext_cdf)]
ext_cdf$SMWT_s <- scale(boxcox_transformation(ext_cdf$SMWT))[1:nrow(ext_cdf)]


## Generate CTG distribution plots

# Un-scaled plots
A1 <- ggplot(data = V2_rec_cdf, aes(x=Mode)) +
  geom_histogram(bins=30, color="black", fill="brown1") +
  scale_x_continuous(breaks=seq(0,1600,400), 
                     limits = c(0,1600)) +
  ggtitle("OPTIMISTIC baseline") +
  xlab("CTG repeat") +
  ylab("Frequency")+
  ptheme

B1 <- ggplot(data = V4_rec_cdf, aes(x=Mode)) +
  geom_histogram(bins=30, color="black", fill="chocolate1") +
  scale_x_continuous(breaks=seq(0,1600,400), 
                     limits = c(0,1600)) +
  ggtitle("OPTIMISTIC 10 months") +
  xlab("CTG repeat") +
  ylab("Frequency")+
  ptheme

C1 <- ggplot(data = ext_cdf, aes(x=CTG)) +
  geom_histogram(bins=30, color="black", fill="cyan2") +
  scale_x_continuous(breaks=seq(0,1600,400), 
                     limits = c(0,1600)) +
  ggtitle("Canadian cohort") +
  xlab("CTG repeat") +
  ylab("Frequency")+
  ptheme

# Scaled plots
A2 <- ggplot(data = V2_rec_cdf, aes(x=Mode_s)) +
  geom_histogram(bins=30, color="black", fill="brown1") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3, 3)) +
  xlab(paste("CTG repeat",
             "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

B2 <- ggplot(data = V4_rec_cdf, aes(x=Mode_s)) +
  geom_histogram(bins=30, color="black", fill="chocolate1") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3, 3)) +
  xlab(paste("CTG repeat",
             "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

C2 <- ggplot(data = ext_cdf, aes(x=CTG_s)) +
  geom_histogram(bins=30, color="black", fill="cyan2") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3, 3)) +
  xlab(paste("CTG repeat",
                 "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

# Arrange CTG plots in grid & save
CTG_plots <- plot_grid(A1, B1, C1, A2, B2, C2,
          ncol=3, nrow=2,
          labels="AUTO",
          label_size = 20,
          align = "hv")
          
ggsave(CTG_plots, 
       file ="SFig1 CTG_distribution_comparisons.png", 
       units="cm", height = 20, width = 30, 
       dpi = 600,
       device = "png")


## Generate SMWT distribution plots

# Un-scaled plots
A1 <- ggplot(data = V2_rec_cdf, aes(x=SMWT)) +
  geom_histogram(bins=30, color="black", fill="brown1") +
  scale_x_continuous(breaks=seq(0,800,200), 
                     limits = c(0,800)) +
  ggtitle("OPTIMISTIC baseline") +
  xlab("6MWT score (m)") +
  ylab("Frequency")+
  ptheme

B1 <- ggplot(data = V4_rec_cdf, aes(x=SMWT)) +
  geom_histogram(bins=30, color="black", fill="chocolate1") +
  scale_x_continuous(breaks=seq(0,800,200), 
                     limits = c(0,800)) +
  ggtitle("OPTIMISTIC 10 months") +
  xlab("6MWT score (m)") +
  ylab("Frequency")+
  ptheme

C1 <- ggplot(data = ext_cdf, aes(x=SMWT)) +
  geom_histogram(bins=30, color="black", fill="cyan2") +
  scale_x_continuous(breaks=seq(0,800,200), 
                     limits = c(0,800)) +
  ggtitle("Canadian cohort") +
  xlab("6MWT score (m)") +
  ylab("Frequency")+
  ptheme

# Scaled plots
A2 <- ggplot(data = V2_rec_cdf, aes(x=SMWT_s)) +
  geom_histogram(bins=30, color="black", fill="brown1") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3.1, 3.1)) +
  xlab(paste("6MWT score",
             "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

B2 <- ggplot(data = V4_rec_cdf, aes(x=SMWT_s)) +
  geom_histogram(bins=30, color="black", fill="chocolate1") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3.1, 3.1)) +
  xlab(paste("6MWT score",
             "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

C2 <- ggplot(data = ext_cdf, aes(x=SMWT_s)) +
  geom_histogram(bins=30, color="black", fill="cyan2") +
  scale_x_continuous(breaks=seq(-3,3,1), 
                     limits = c(-3.1, 3.1)) +
  xlab(paste("6MWT score",
             "[Box-Cox transformed z-score]", sep="\n")) +
  ylab("Frequency")+
  ptheme

# Arrange CTG plots in grid & save
SMWT_plots <- plot_grid(A1, B1, C1, A2, B2, C2,
                       ncol=3, nrow=2,
                       labels="AUTO",
                       label_size = 20,
                       align = "hv")

ggsave(SMWT_plots, 
       file ="SFig2 SMWT_distribution_comparisons.png", 
       units="cm", height = 20, width = 30, 
       dpi = 600,
       device = "png")




