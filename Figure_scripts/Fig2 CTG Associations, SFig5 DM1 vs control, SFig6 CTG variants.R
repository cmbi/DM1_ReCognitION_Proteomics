
## Script generated on 29/7/24, last update 02/06/25
## Figure 2: Top CTG-repeat associations, ReCognitION + External validation; 
## SFigure 5: Top CTG and SMWT associations vs controls
## SFigure 6: Impact of CTG-repeat interruptions on CTG-associations


####################
## Load libraries ##
####################
library(ggplot2)       #V 3.3.6
library(cowplot)       #V 1.1.1
library(grid)          #V 4.0.4
library(ggpubr)        #V 0.4.0


################
## Load files ##
################

## External validation proteomics data
load(file = "all_p_ext_val_res.RDATA")
cdf_ext <- cdf
rm(c, cdf, CTG_res_df, SMWT_res_df)

## ReCognitION proteomics statistical results
load (file = "DIA-NN_p_Protein_statistics.RDATA")
cdf_rec <- cdf 
rm(cdf)

## Table 2 to obtain overlapping significant/validated protein groups
load(file = "Table2.RDATA")

## External proteomics data, normalized and log2 scaled, no minimum threshold
load(file = "ext_val_normalized_proteins_no_minimum.RData")

## Set ggplot2 theme
ptheme <- theme(
  aspect.ratio = 1,
  plot.title = element_text(color = "black", size = 16, face="bold"),
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


#################################################
## Figure 2: Protein - CTG repeat associations ##
#################################################

### Volcano plot of CTG effect
## Subset relevant data for plotting
df <- data.frame(protein = rownames(Mode_res_df),
                 p.value = Mode_res_df$`Mode_Pr(>|t|)`,
                 FDR = Mode_res_df$Mode_FDR,
                 effect = Mode_res_df$Mode_Estimate,
                 pcor = Mode_res_df$pcor)

## Add significance labels, including for external validation, order for plotting
sum(is.na(df$FDR))
df$sig <- ifelse(df$FDR < 0.05, "sig", "not_sig")
for (protein in df$protein){
  if (protein %in% T2$Rec_PG[T2$OM=="CTG"]){
    df$sig[df$protein == protein] <- "sig_ext"
  }
}
df <- df[order(df$sig),]

## Generate plot
A <- ggplot(df, aes(x = effect*100, y = -log10(p.value), color=sig))+
  geom_point(size= 1) +
  scale_color_manual(breaks = c("not_sig", "sig", "sig_ext"),
                     values = c("grey", "black", "red")) +
  xlab(paste("CTG repeat effect size",
             "(per 100 CTGs)", sep="\n")) +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 22), 
                     breaks = seq(0, 20, 4),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-0.15, 0.15), 
                     breaks = seq(-0.15, 0.15, 0.15)) +
  ptheme +
  theme(legend.position = "none",
        legend.box.just = "right", 
        aspect.ratio = 2)


### Example plots: top 2 validated hits
## Obtain top 2 hits based on absolute Pearson Rho
df_ex <- df[df$sig == "sig_ext",]
df_ex <- df_ex[order(abs(df_ex$pcor), decreasing=T),]
df_ex$protein[1:2]  # P01857;P01860;P01861 = IGHG1;IGHG3;IGHG4 | P01024 = C3

## OPTIMISTIC/ReCognitION plots
# P01857;P01860;P01861 = IGHG1;IGHG3;IGHG4
pcor_base <- round(cor(cdf_rec[cdf_rec$Visit=="V2", "Mode"], 
                 cdf_rec[cdf_rec$Visit=="V2", "P01857;P01860;P01861"], 
                 use="pairwise.complete.obs"),2)
pcor_10M <- round(cor(cdf_rec[cdf_rec$Visit=="V4", "Mode"], 
                      cdf_rec[cdf_rec$Visit=="V4", "P01857;P01860;P01861"], 
                      use="pairwise.complete.obs"),2)

B <- ggplot(cdf_rec, aes(x=Mode, y=`P01857;P01860;P01861`, color=Visit))+
  scale_color_manual(breaks = c("V2", "V4"),
                     values = c("red", "blue"),
                     labels = c("V2" = paste0("Baseline, Rho = ", pcor_base),
                                "V4" = paste0("10 Months, Rho = ", pcor_10M)))+
  ggtitle("OPTIMISTIC IGHG1/IGHG3/IGHG4") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("PG Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(-1.4, 5)) +
  ptheme +
  theme(
    legend.position = c(0.60, 0.88),
    legend.box.just = "right",
    legend.title = element_blank()
  )

# P01024 = C3
pcor_base <- round(cor(cdf_rec[cdf_rec$Visit=="V2", "Mode"], 
                       cdf_rec[cdf_rec$Visit=="V2", "P01024"], 
                       use="pairwise.complete.obs"),2)
pcor_10M <- round(cor(cdf_rec[cdf_rec$Visit=="V4", "Mode"], 
                      cdf_rec[cdf_rec$Visit=="V4", "P01024"], 
                      use="pairwise.complete.obs"),2)

C <- ggplot(cdf_rec, aes(x=Mode, y=P01024, color=Visit))+
  scale_color_manual(breaks = c("V2", "V4"),
                     values = c("red", "blue"),
                     labels = c("V2" = paste0("Baseline, Rho = ", pcor_base),
                                "V4" = paste0("10 Months, Rho = ", pcor_10M)))+
  ggtitle("OPTIMISTIC C3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(3.5, 6.5)) +
  ptheme +
  theme(
    legend.position = c(0.6, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
    )


## External validation plots
# P01860 = IGHG3
pcor_ext <- round(cor(cdf_ext$CTG, cdf_ext$P01860, use="pairwise.complete.obs"),2)

D <- ggplot(cdf_ext, aes(x=CTG, y=P01860))+
  ggtitle("Validation IGHG3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1400,400), 
                     limits = c(0,1600)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.6, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )

# P01024 = C3
pcor_ext <- round(cor(cdf_ext$CTG, cdf_ext$P01024, use="pairwise.complete.obs"),2)

E <- ggplot(cdf_ext, aes(x=CTG, y=P01024))+
  ggtitle("Validation C3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("CTG repeat") +
  ylab("PG Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1400,400), 
                     limits = c(0,1600)) +
  scale_y_continuous(limits = c(3.8, 5.5)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.6, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )

## Arrange plots and save
p1 <- plot_grid(A,
                labels=c("A"),
                label_size = 20)

p2 <- plot_grid(B, D, C, E,
                ncol = 2,
                align ="hv",
                labels = c("B", "C", "D", "E"),
                label_size = 20)

Fig2 <- plot_grid(p1, p2, nrow=1,
                  rel_widths = c(1, 2))

ggsave(Fig2, 
       file ="Fig2_CTG_associations.png", 
       units="cm", height = 20, width = 35, 
       dpi = 600,
       device = "png")


#####################################################
## SFigure 5: Top CTG and 6MWT hits DM1 vs control ##
#####################################################


## Transverse dataframe and label cases
all_p_avg_t <- as.data.frame(t(all_p_avg))
all_p_avg_t$Cohort <- c(rep("Canadian DM1", 56), 
                         rep("German Control", 10),
                         rep("German DM1", 13))

all_p_avg_t$Patient <- c(rep("DM1", 56), 
                         rep("Control", 10),
                         rep("DM1", 13))
  
## SFigure 6: first row = Top 4 CTG-immunoglobulin hits DM1 vs control
# Needs dataframe generated above (Figure 2)
T2[T2$OM=="CTG",][1:5,c("Rec_PG", "Ext_PG")]


# P01857;P01860;P01861 - matched with: P01860 / IGHG3
A <- ggplot(all_p_avg_t, aes(x=Patient, y=P01860))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("IGHG3") +
  xlab("") +
  scale_y_continuous(limits = c(2, 7)) +
  ylab("Protein Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 6.85, label.x = 0.9)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# B9A064 = IGLL5, matched with B9A064;P0CG04;P0DOX8;P0DOY2;P0DOY3 = IGLL5/IGLC(1,2,3)/IGL1        
title_B <- paste("IGLC(1,2,3)", "IGL1/IGLL5", sep="\n")
B <- ggplot(all_p_avg_t, aes(x=Patient, y=`B9A064;P0CG04;P0DOX8;P0DOY2;P0DOY3`))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle(title_B) +
  xlab("") +
  scale_y_continuous(limits = c(2, 5)) +
  ylab("PG Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 4.9)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# A0M8Q6;B9A064 = IGLC7;IGLL5 matched with A0M8Q6;B9A064;P0CF74;P0CG04;P0DOX8;P0DOY2;P0DOY3 = IGLL5/IGL1/IGLC(1,2,3,6,7)
C <- ggplot(all_p_avg_t, aes(x=Patient, y=`A0M8Q6;B9A064;P0CF74;P0CG04;P0DOX8;P0DOY2;P0DOY3`))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle(title_c) +
  xlab("") +
  ylab("PG Intensity (log2)") +
  scale_y_continuous(limits = c(2, 5.2)) +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 5.1)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# A0A075B6K5;P80748 = IGLV3-9/IGLV3-21 exact matched with A0A075B6K5;P80748
D <- ggplot(all_p_avg_t, aes(x=Patient, y=`A0A075B6K5;P80748`))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("IGLV3-9/IGLV3-21") +
  xlab("") +
  scale_y_continuous(limits = c(-2.5, 1)) +
  ylab("PG Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 0.9)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))


## Second row: top 3 SMWT - complement hits and SMWT-ITIH3 
T2[T2$OM=="6MWT",][1:5,c("Rec_PG", "Ext_PG")]

#P01024 = C3
E <- ggplot(all_p_avg_t, aes(x=Patient, y=P01024))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("C3") +
  xlab("") +
  scale_y_continuous(limits = c(5.3, 7.5)) +
  ylab("Protein Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 7.4)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# P05156 = CFI
plotF <- ggplot(all_p_avg_t, aes(x=Patient, y=P05156))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("CFI") +
  xlab("") +
  scale_y_continuous(limits = c(1.4, 3.3)) +
  ylab("Protein Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 3.24)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# P01031 = C5
G <- ggplot(all_p_avg_t, aes(x=Patient, y=P01031))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("C5") +
  xlab("") +
  scale_y_continuous(limits = c(2, 3.6)) +
  ylab("Protein Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 3.54)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

# Q06033 = ITIH3
H <- ggplot(all_p_avg_t, aes(x=Patient, y=Q06033))+
  geom_boxplot(fill='#A4A4A4')+
  ggtitle("ITIH3") +
  xlab("") +
  scale_y_continuous(limits = c(0.2, 2.9)) +
  ylab("Protein Intensity (log2)") +
  stat_compare_means(method ="t.test", paired = F, size=5, label.y = 2.8)+
  ptheme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(face="bold", color = "black", size = 16))

## Arrange and save
SFig5 <- plot_grid(A, B, C, D, 
                   E, plotF, G, H,
                   nrow = 2,
                   align ="hv",
                   labels = c("AUTO"),
                   label_size = 20)

ggsave(SFig5, 
       file ="SFig5_DM1_vs_Control.png", 
       units="cm", height = 20, width = 35, 
       dpi = 600,
       device = "png")


###################################################
## SFigure 6: Effect of CTG-repeat interruptions ##
###################################################

## Needs dataframe generated above (Figure 2)
# P01857;P01860;P01861 = IGHG1/IGHG3/IGHG4 
# B9A064 = IGLL5
# A0M8Q6;B9A064 = IGLC7;IGLL5
# P01024 = C3
# P05156 = CFI
# P08603 = CFH


## Subset cases where VariantRepeats have been estimated
cdf_rec_var <- cdf_rec[!cdf_rec$VariantRepeats=="NA",]
table(cdf_rec_var$VariantRepeats)
cdf_rec_var$VariantRepeats <- ifelse(cdf_rec_var$VariantRepeats==0, "no", "yes")

## P01857;P01860;P01861 = IGHG1/IGHG3/IGHG4 
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "P01857;P01860;P01861"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "P01857;P01860;P01861"], 
                     use="pairwise.complete.obs"),2)

A <- ggplot(cdf_rec_var, aes(x=Mode, y=`P01857;P01860;P01861`, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("IGHG1/IGHG3/IGHG4") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("PG Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(-1.4, 6)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

## B9A064 = IGLL5
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "B9A064"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "B9A064"], 
                     use="pairwise.complete.obs"),2)

B <- ggplot(cdf_rec_var, aes(x=Mode, y=B9A064, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("IGLL5") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(2, 7)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )


## A0M8Q6;B9A064 = IGLC7/IGLL5
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "A0M8Q6;B9A064"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "A0M8Q6;B9A064"], 
                     use="pairwise.complete.obs"),2)

C <- ggplot(cdf_rec_var, aes(x=Mode, y=`A0M8Q6;B9A064`, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("IGLC7/IGLL5") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("PG Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(-0.6, 6)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )


## P01024 = C03
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "P01024"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "P01024"], 
                     use="pairwise.complete.obs"),2)

D <- ggplot(cdf_rec_var, aes(x=Mode, y=P01024, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("C3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(3.5, 6)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

## P05156 = CFI
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "P05156"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "P05156"], 
                     use="pairwise.complete.obs"),2)

E <- ggplot(cdf_rec_var, aes(x=Mode, y=P05156, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("CFI") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(-1.5, 2.5)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

## P08603 = CFH
pcor_noVR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "Mode"], 
                       cdf_rec_var[cdf_rec_var$VariantRepeats=="no", "P08603"], 
                       use="pairwise.complete.obs"),2)
pcor_VR <- round(cor(cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "Mode"], 
                     cdf_rec_var[cdf_rec_var$VariantRepeats=="yes", "P08603"], 
                     use="pairwise.complete.obs"),2)

plotF <- ggplot(cdf_rec_var, aes(x=Mode, y=P08603, color=VariantRepeats))+
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("red", "blue"),
                     labels = c("no" = paste0("No Interruption, Rho = ", pcor_noVR),
                                "yes" = paste0("CTG Interruption, Rho = ", pcor_VR)))+
  ggtitle("CFH") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1200,400), 
                     limits = c(0,1400)) +
  scale_y_continuous(limits = c(1.5, 5)) +
  ptheme +
  theme(
    legend.position = c(0.5, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

## Arrange and safe
SFig6 <- plot_grid(A, B, C, D, E, plotF,
                   ncol = 3,
                   align ="hv",
                   labels = c("AUTO"),
                   label_size = 20)

ggsave(SFig6, 
       file ="SFig6_CTG_Interruptions.png", 
       units="cm", height = 20, width = 35, 
       dpi = 600,
       device = "png")


