
## Script generated on 29/7/24, last update 12/05/25
## Figure 3: Top 6MWT-repeat associations, ReCognitION + External validation; 
## Figure 4: ELISA based validation of ITIH3

####################
## Load libraries ##
####################
library(ggplot2)       #V 3.3.6
library(cowplot)       #V 1.1.1
library(grid)          #V 4.0.4


################
## Load files ##
################

## External validation proteomics statistical results
load(file = "all_p_ext_val_res.RDATA")
cdf_ext <- cdf
rm(cdf, CTG_res_df, SMWT_res_df)

## ReCognitION proteomics statistical results
load (file = "DIA-NN_p_Protein_statistics.RDATA")
cdf_rec <- cdf 
rm(cdf)

## Table 2 to obtain overlapping significant/validated protein groups
load(file = "Table2.RDATA")

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
## Figure 3: Protein - 6MWT repeat associations ##
#################################################

### Volcano plot of 6MWT effect
## Subset relevant data for plotting
SMWT_rec <- outcome_coef[["SMWT"]]

df <- data.frame(protein = rownames(SMWT_rec),
                 p.value = SMWT_rec$`SMWT_Pr(>|t|)`,
                 FDR = SMWT_rec$SMWT_FDR,
                 effect = SMWT_rec$SMWT_Estimate,
                 pcor = SMWT_rec$pcor)

## Add significance labels, including for external validation, order for plotting
sum(is.na(df$FDR))
df$sig <- ifelse(df$FDR < 0.05, "sig", "not_sig")
for (protein in df$protein){
  if (protein %in% T2$Rec_PG[T2$OM=="6MWT"]){
    df$sig[df$protein == protein] <- "sig_ext"
  }
}
df <- df[order(df$sig),]

## Generate plot
A <- ggplot(df, aes(x = effect*100, y = -log10(p.value), color=sig))+
  geom_point(size= 1) +
  scale_color_manual(breaks = c("not_sig", "sig", "sig_ext"),
                     values = c("grey", "black", "red")) +
  xlab(paste("6MWT score effect size",
             "(per 100m)", sep="\n")) +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 10), 
                     breaks = seq(0, 10, 2.5),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-0.3, 0.3), 
                     breaks = seq(-0.3, 0.3, 0.3)) +
  ptheme +
  theme(legend.position = "none",
        legend.box.just = "right", 
        aspect.ratio = 2)


### Example plots: top 2 validated hits
## Obtain top 2 hits
df_ex <- df[df$sig == "sig_ext",]
df_ex <- df_ex[order(abs(df_ex$pcor), decreasing=T),]
df_ex$protein[1:2]  # P01024= C3; P05156 = CFI.

## OPTIMISTIC/ReCognitION plots
# P01024= C3
pcor_base <- round(cor(cdf_rec[cdf_rec$Visit=="V2", "SMWT"], 
                 cdf_rec[cdf_rec$Visit=="V2", "P01024"], 
                 use="pairwise.complete.obs"),2)
pcor_10M <- round(cor(cdf_rec[cdf_rec$Visit=="V4", "SMWT"], 
                      cdf_rec[cdf_rec$Visit=="V4", "P01024"], 
                      use="pairwise.complete.obs"),2)

B <- ggplot(cdf_rec, aes(x=SMWT, y=P01024, color=Visit))+
  scale_color_manual(breaks = c("V2", "V4"),
                     values = c("red", "blue"),
                     labels = c("V2" = paste0("Baseline, Rho = ", pcor_base),
                                "V4" = paste0("10 Months, Rho = ", pcor_10M)))+
  ggtitle("OPTIMISTIC C3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  scale_y_continuous(limits = c(3.5, 6)) +
  ptheme +
  theme(
    legend.position = c(0.60, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

# P05156 = CFI
pcor_base <- round(cor(cdf_rec[cdf_rec$Visit=="V2", "SMWT"], 
                       cdf_rec[cdf_rec$Visit=="V2", "P05156"], 
                       use="pairwise.complete.obs"),2)
pcor_10M <- round(cor(cdf_rec[cdf_rec$Visit=="V4", "SMWT"], 
                      cdf_rec[cdf_rec$Visit=="V4", "P05156"], 
                      use="pairwise.complete.obs"),2)

C <- ggplot(cdf_rec, aes(x=SMWT, y=`P05156`, color=Visit))+
  scale_color_manual(breaks = c("V2", "V4"),
                     values = c("red", "blue"),
                     labels = c("V2" = paste0("Baseline, Rho = ", pcor_base),
                                "V4" = paste0("10 Months, Rho = ", pcor_10M)))+
  ggtitle("OPTIMISTIC CFI") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  scale_y_continuous(limits = c(-1.3, 2.2)) +
  ptheme +
  theme(
    legend.position = c(0.60, 0.85),
    legend.box.just = "right",
    legend.title = element_blank()
  )

## External validation plots
# P01024= C3;
pcor_ext <- round(cor(cdf_ext$SMWT, cdf_ext$P01024, use="pairwise.complete.obs"),2)

D <- ggplot(cdf_ext, aes(x=SMWT, y=P01024))+
  ggtitle("Validation C3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.6, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )


# P05156 = CFI
pcor_ext <- round(cor(cdf_ext$SMWT, cdf_ext$P05156, use="pairwise.complete.obs"),2)

E <- ggplot(cdf_ext, aes(x=SMWT, y=P05156))+
  ggtitle("Validation CFI") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
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

Fig3 <- plot_grid(p1, p2, nrow=1,
                  rel_widths = c(1, 1.8))

ggsave(Fig3, 
       file ="Fig3_SMWT_associations.png", 
       units="cm", height = 20, width = 35, 
       dpi = 600,
       device = "png")


#############################################
## Figure 4 (ELISA based) ITIH3 validation ##
#############################################

## Panel A: ITIH3 vs 6MWT ReCognitION [baseline vs 10 months]
pcor_base <- round(cor(cdf_rec[cdf_rec$Visit=="V2", "SMWT"], 
                       cdf_rec[cdf_rec$Visit=="V2", "Q06033"], 
                       use="pairwise.complete.obs"),2)
pcor_10M <- round(cor(cdf_rec[cdf_rec$Visit=="V4", "SMWT"], 
                      cdf_rec[cdf_rec$Visit=="V4", "Q06033"], 
                      use="pairwise.complete.obs"),2)

A <- ggplot(cdf_rec, aes(x=SMWT, y=Q06033, color=Visit))+
  scale_color_manual(breaks = c("V2", "V4"),
                     values = c("red", "blue"),
                     labels = c("V2" = paste0("Baseline, Rho = ", pcor_base),
                                "V4" = paste0("10 Months, Rho = ", pcor_10M)))+
  ggtitle("OPTIMISTIC ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  ptheme +
  theme(
    legend.position = c(0.60, 0.88),
    legend.box.just = "right",
    legend.title = element_blank(),
  )

## Panel B: ITIHI3 vs 6MWT External validation
pcor_ext <- round(cor(cdf_ext$SMWT, cdf_ext$Q06033, use="pairwise.complete.obs"),2)

B <- ggplot(cdf_ext, aes(x=SMWT, y=Q06033))+
  ggtitle("Validation ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("6MWT score (m)") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  scale_y_continuous(limits = c(-1.4, 1.2)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.6, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )

## Panel C: ELISA based ITIH3 vs 6MWT External validation 
pcor_ext <- round(cor(c$SMWT, c$ITIH3log2, use="pairwise.complete.obs"),2)

C <- ggplot(c, aes(x=SMWT, y=ITIH3log2))+
  ggtitle("ELISA Based Validation ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("6MWT score (m)") +
  ylab("ITIH3 (log2 mg/ml)") +
  scale_x_continuous(breaks=seq(0,750,250), 
                     limits = c(0,850)) +
  scale_y_continuous(limits = c(-5, 1)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.6, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )

## Panel D: ELISA based ITIH3 vs Grip strength (left hand) External validation
pcor_ext <- round(cor(c$ITIH3log2, c$gripG, use="pairwise.complete.obs"),2)

D <- ggplot(c, aes(x=gripG, y=ITIH3log2))+
  ggtitle("ELISA Based Validation ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("Grip strengh left hand (kg)") +
  ylab("ITIH3 (log2 mg/ml)") +
  scale_x_continuous(breaks=seq(0,50,12.5), 
                     limits = c(0,50)) +
  scale_y_continuous(limits = c(-5, 1)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.58, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )


## Panel E: ELISA based ITIH3 vs Grip strength (right hand) External validation 
pcor_ext <- round(cor(c$ITIH3log2, c$gripD, use="pairwise.complete.obs"),2)

E <- ggplot(c, aes(x=gripD, y=ITIH3log2))+
  ggtitle("ELISA Based Validation ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col="black") +
  xlab("Grip strengh right hand (kg)") +
  ylab("ITIH3 (log2 mg/ml)") +
  scale_x_continuous(breaks=seq(0,50,12.5), 
                     limits = c(0,50)) +
  scale_y_continuous(limits = c(-5, 1)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor_ext, 2)), 
    x=0.58, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  ptheme +
  theme(
    legend.box.just = "right",
    legend.title = element_blank()
  )


## Panel F: ITIH3 vs CTG-repeat in both studies
# Subset relevant data from both studies & merge
cdf_rec_ITIH3 <- cdf_rec[, c("Mode", "Q06033")]
colnames(cdf_rec_ITIH3)[1] <- "CTG"
cdf_rec_ITIH3$label <- rep("OPTIMISTIC", nrow(cdf_rec_ITIH3))
rec_ITIH3_CTG_cor <- round(cor(cdf_rec_ITIH3$CTG, cdf_rec_ITIH3$Q06033, 
                         use="pairwise.complete.obs"),2)

cdf_ext_ITIH3 <- cdf_ext[, c("CTG", "Q06033")]
cdf_ext_ITIH3$label <- rep("Validation", nrow(cdf_ext_ITIH3))
ext_ITIH3_CTG_cor <- round(cor(cdf_ext_ITIH3$CTG, cdf_ext_ITIH3$Q06033, 
                         use="pairwise.complete.obs"),2)

table(colnames(cdf_rec_ITIH3) == colnames(cdf_ext_ITIH3))
fdf <- as.data.frame(rbind(cdf_rec_ITIH3, cdf_ext_ITIH3))

plotF <- ggplot(fdf, aes(x=CTG, y=Q06033, color=label))+
  scale_color_manual(breaks = c("OPTIMISTIC", "Validation"),
                     values = c("red", "blue"),
                     labels = c("OPTIMISTIC" = paste0("OPTIMISTIC, Rho = ", 
                                                      rec_ITIH3_CTG_cor),
                                "Validation" = paste0("Validation, Rho = ", 
                                                      ext_ITIH3_CTG_cor)))+
  ggtitle("Both Cohorts ITIH3") +
  geom_point(size=1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F) +
  xlab("CTG repeat") +
  ylab("Protein Intensity (log2)") +
  scale_x_continuous(breaks=seq(0,1600,400), 
                     limits = c(0,1600)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  ptheme +
  theme(
    legend.position = c(0.60, 0.86),
    legend.box.just = "right",
    legend.title = element_blank()
  )


## Arrange plots and save
Fig4 <- plot_grid(A, B, C, D, E, plotF,
                ncol = 3,
                align ="hv",
                labels = "AUTO",
                label_size = 20)

ggsave(Fig4, 
       file ="Fig4_ITIH3.png", 
       units="cm", height = 20, width = 35, 
       dpi = 600,
       device = "png")


