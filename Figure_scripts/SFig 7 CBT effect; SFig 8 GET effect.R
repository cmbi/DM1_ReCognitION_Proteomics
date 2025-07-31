
## Script generated on 29/7/24, last update 25/04/25
## SFigure 7: CBT effect
## SFigure 8: GET effect

####################
## Load libraries ##
####################
library(ggplot2)       #V 3.3.6
library(cowplot)       #V 1.1.1
library(grid)          #V 4.0.4


################
## Load files ##
################

## ReCognitION proteomics statistical results
load (file = "DIA-NN_p_Protein_statistics.RDATA")

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


###########################
## SFigure 7: CBT effect ## 
###########################

## Volcano plot of GET effect
df <- data.frame(protein = rownames(Vis_Treat_res_df),
                 p.value = Vis_Treat_res_df$`Visit_Treatment_Pr(>|t|)`,
                 FDR = Vis_Treat_res_df$Visit_Treatment_FDR,
                 effect = Vis_Treat_res_df$Visit_Treatment_Estimate)

# Add significance labels, order for plotting
sum(is.na(df$FDR))
df$sig <- ifelse(df$FDR < 0.05, "sig", "not_sig")
df <- df[order(df$sig),]

A <- ggplot(df, aes(x = effect, y = -log10(p.value), color=sig))+
  geom_point(size= 1) +
  scale_color_manual(breaks = c("not_sig", "sig"),
                     values = c("grey", "black")) +
  xlab("CBT effect size") +
  ylab("-10log(p-value)") +
  scale_x_continuous(limits = c(-0.5, 0.5), 
                     breaks = seq(-0.5, 0.5, 0.5)) +
  ptheme +
  theme(aspect.ratio = 2,
        legend.position = "none")

## Example plots
cdf$group <- cdf$TreatmentCode
sum(is.na(cdf$group))
cdf$group <- ifelse(cdf$group == 0, "Control", "CBT")

cdf$Timepoint <- cdf$Visit
sum(is.na(cdf$Timepoint))
cdf$Timepoint <- ifelse(cdf$Timepoint == "V4", "Follow-up 10M", "Baseline")

Vis_Treat_res_df <- Vis_Treat_res_df[order(Vis_Treat_res_df$`Visit_Treatment_Pr(>|t|)`),]
table(Vis_Treat_res_df$Visit_Treatment_FDR < 0.05)
rownames(Vis_Treat_res_df)[1:8]
CBT_names <- c("GC", "C5", "CP", "IGHG3", "C9", "IGHM", "ITIH4", "CPB2")

CBT_box_list <- list()
for (i in 1:8){
  pname <- rownames(Vis_Treat_res_df)[i]
  CBT_box_list[[pname]] <- ggplot(cdf, aes_string(y=cdf[,pname], x="group", fill="Timepoint"))+
    ggtitle(CBT_names[i])+
    geom_boxplot() +
    scale_fill_manual(values=c('#999999', '#E69F00')) +
    xlab("") +
    ylab("") +
    ptheme +
    theme(legend.position = "none",
          axis.text.x = element_text(face="bold", color = "black", size = 16)
    )
}

# Generate one plot separately with a legend and store legend in list
plot_ex <- ggplot(cdf, aes_string(y=cdf[,pname], x="group", fill="Timepoint"))+
  ggtitle(CBT_names[5])+
  geom_boxplot() +
  scale_fill_manual(values=c('#999999', '#E69F00')) +
  xlab("") +
  ylab("") +
  ptheme +
  theme(
    legend.justification = "top",
    legend.title = element_text(size=20, face="bold"),
    legend.text = element_text(face="bold",color = "black", size = 20),)

CBT_box_list[[9]] <- get_legend(plot_ex)


## Cast subplots into one and save
P1 <- plot_grid(A, labels = c("A"),
                label_size = 20)


P2 <- plot_grid(plotlist = CBT_box_list, ncol=3,
                labels= c("B", "C", "D", "E", "F", "G", "H", "I"),
                label_size = 20,
                align ="h",
                axis="t")


SFig7 <- plot_grid(P1, P2, nrow=1,
                   rel_widths = c(0.5,1))

ggsave(SFig7, 
       file ="Figures/SFig7_CBT_effect.png", 
       units="cm", height = 25, width = 40, 
       dpi = 600,
       device = "png")


###########################
## SFigure 8: GET effect ##
###########################

## Volcano plot of GET effect
df <- data.frame(protein = rownames(Vis_Get_res_df),
                 p.value = Vis_Get_res_df$`Visit_Get_Pr(>|t|)`,
                 FDR = Vis_Get_res_df$Visit_Get_FDR,
                 effect = Vis_Get_res_df$Visit_Get_Estimate)

# Add significance labels, order for plotting
sum(is.na(df$FDR))
df$sig <- ifelse(df$FDR < 0.05, "sig", "not_sig")
df <- df[order(df$sig),]

A <- ggplot(df, aes(x = effect, y = -log10(p.value), color=sig))+
  geom_point(size= 1) +
  scale_color_manual(breaks = c("not_sig", "sig"),
                     values = c("grey", "black")) +
  xlab("GET effect size") +
  ylab("-10log(p-value)") +
  scale_x_continuous(limits = c(-1.3, 1.3), 
                     breaks = seq(-1, 1, 1)) +
  ptheme +
  theme(aspect.ratio = 2,
        legend.position = "none")

## Example plots
cdf$group <- cdf$GradedExerciseTherapy
sum(is.na(cdf$group))
cdf$group <- ifelse(cdf$group == 0, "Control", "GET")

cdf$Timepoint <- cdf$Visit
sum(is.na(cdf$Timepoint))
cdf$Timepoint <- ifelse(cdf$Timepoint == "V4", "Follow-up 10M", "Baseline")

Vis_Get_res_df <- Vis_Get_res_df[order(Vis_Get_res_df$`Visit_Get_Pr(>|t|)`),]
sum(Vis_Get_res_df$Visit_Get_FDR < 0.05)
rownames(Vis_Get_res_df)[1:2] #P33151 = CADH5; P02751 = FN1

#P33151 = Cadherin-5
B <- ggplot(cdf, aes_string(y=cdf[,"P33151"], x="group", fill="Timepoint"))+
  ggtitle("Cadherin-5")+
  geom_boxplot() +
  scale_fill_manual(values=c('#999999', '#E69F00')) +
  xlab("") +
  ylab("") +
  ptheme +
  theme(legend.justification = "top",
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(face="bold",color = "black", size = 20),
        axis.text.x = element_text(face="bold", color = "black", size = 16))

## P02751 = Fibronectin
C <- ggplot(cdf, aes_string(y=cdf[,"P02751"], x="group", fill="Timepoint"))+
  ggtitle("Fibronectin")+
  geom_boxplot() +
  scale_fill_manual(values=c('#999999', '#E69F00')) +
  xlab("") +
  ylab("") +
  ptheme +
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold", color = "black", size = 16))


## Cast subplots into one and save
P1 <- plot_grid(A, labels = c("A"),
                label_size = 20)
P2 <- plot_grid(B, C, 
                ncol=1,
                labels= c("B", "C"),
                label_size = 20,
                align = "v",
                axis = "r")
SFig8 <- plot_grid(P1, P2, nrow=1,
                   rel_widths = c(1,1.5))

ggsave(SFig8, 
       file ="SFig8_GET_effect.png", 
       units="cm", height = 20, width = 25, 
       dpi = 600,
       device = "png")
