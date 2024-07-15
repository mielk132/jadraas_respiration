#Jädraås Respiration Manuscript


#Author: Louis Mielke
#Date Updated: 2021 MAY 28

#Supplementary Figure 3 

# (a) Proportion of ericaceous shrub DNA out of total ITS2 reads 
# (b) and proportion of the ITS2 gene region from ectomycorrhizal fungi out of total fungal ITS2 reads
# data from 50 µm-mesh bags with oven dried humus from the same site
# incubated for 17 months from June 2018 to November 2019.
# The percent relative abundance (%) from Ericaceae roots and EcM mycelium are calculated as total number of reads for those organisms out of the total ITS2 reads and total fungal ITS2 reads, respectively. 
# Means ± standard errors (n=8), besides (n=7) for the pine root removal. 
# Significant differences among treatments are indicated with different letters (P < 0.05).


rm(list=ls()) #clear the workspace


#rabund from all_clusters_scata4703_25MAY2020.xlsx file
library(tidyverse)
library(nlme)
library(ggpattern)
library(ggplot2)
library(gridExtra)

setwd("/Users/lske0002/Projects/mycorrhizal removal/MeshBags/data/humus mass loss/data/set b incubation 17 months")
rabund <- read.csv("Jädraås_Humus_green_2018-2019_rabund.csv")
rabund$block  <- as.factor(rabund$block)
str(rabund)
rabund2 <- na.omit(rabund) %>%
  group_by(treatment) %>%
  drop_na() %>% 
  summarize_at(c("EcM_relative_abund","shrub_relative_abund"), funs(mean, min, max, sd, n()))

#standard error calculation
rabund2$EcM_rabund_se <- rabund2$EcM_relative_abund_sd / sqrt(rabund2$EcM_relative_abund_n)
rabund2$shrub_rabund_se <- rabund2$shrub_relative_abund_sd / sqrt(rabund2$shrub_relative_abund_n)

#ericaceous shrub and plant DNA (2.8% of reads) removed before assessing relative abundance of ecm

ecm <- ggplot(data = rabund2) +
  aes(x = treatment, y = EcM_relative_abund_mean, fill = treatment) +
  geom_bar(stat="identity", color="black", position=position_dodge(),show.legend = FALSE) +
  geom_errorbar(aes(ymin=EcM_relative_abund_mean-EcM_rabund_se, ymax=EcM_relative_abund_mean+EcM_rabund_se), width=.2, position=position_dodge(.9))+
  theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.position = "none",
        axis.line = element_line( size = 1, linetype = "solid"), 
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(colour = "black",size = 12),
        axis.text.y = element_text(colour = "black",size = 12))+
  scale_fill_manual(values = c("#E69F00","#D55E00","#56B4E9", "#009E73","#F0E442"))+
  ylab("ectomycorrhizal relative abundance")+
  xlab("treatment")+
  scale_y_continuous(limits = c(0,.5),expand = (c(0,0)))+
  scale_x_discrete(limits = c("C","DC","E","T","TE"),labels = c("C","DC","SR","PE","SR + PE"))
ecm 

ericaceae <- ggplot(data = rabund2) +
  aes(x = treatment, y = shrub_relative_abund_mean, fill = treatment) +
  geom_bar(stat="identity", color="black", position=position_dodge(), show.legend = FALSE) +
  geom_errorbar(aes(ymin=shrub_relative_abund_mean-shrub_rabund_se, ymax=shrub_relative_abund_mean+shrub_rabund_se), width=.2, position=position_dodge(.9))+
  theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        legend.position = "none",
        axis.line = element_line( size = 1, linetype = "solid"), 
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(colour = "black",size = 12),
        axis.text.y = element_text(colour = "black",size = 12))+
  scale_fill_manual(values = c("#E69F00","#D55E00","#56B4E9", "#009E73","#F0E442"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.10))+
  scale_x_discrete(limits = c("C","DC","E","T","TE"),labels = c("C","DC","SR","PE","SR + PE"))+
  ylab("ericaceae relative abundance")+
  xlab("treatment")
ericaceae

#stats

lme.ecm <- lme(sqrt(EcM_relative_abund) ~ treatment, random = ~ 1|block, data = rabund) 
anova(lme.ecm)
plot(resid(lme.ecm))
emmeans.ecm <- emmeans(lme.ecm, c("treatment"), type = "response")
pairs(emmeans.ecm, type = "response", adjust = "bonf")

lme.eri <- lme(sqrt(shrub_relative_abund) ~ treatment, random = ~ 1|block, data = rabund) 
anova(lme.eri)
plot(resid(lme.eri))
emmeans.eri <- emmeans(lme.eri, c("treatment"), type = "response")
pairs(emmeans.eri, type = "response", adjust = "bonf")

#combine the two rabund plots together
require(gridExtra)
plot1 <- ericaceae
plot2 <- ecm
grid.arrange(plot1, plot2, ncol=2)

setwd("/Users/lske0002/Projects/mycorrhizal removal/Writing/Respiration Paper/Final Graphs")

# Open a pdf file
pdf("rplot_Eri_EcM.pdf") 

grid.arrange(ericaceae, ecm, ncol=2)

# Close the pdf file
dev.off()
