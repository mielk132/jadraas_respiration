#Jädraås Respiration Study
#Table 2
#Soil Respiration and Temp & Moisture Values

#Author: Louis Mielke
#Date Created: 2020 May 28
#Date Updated: 2021 Oct 28


rm(list=ls()) #clear the workspace

library(readxl)
library(ggplot2)
library(ggsci)
library(tidyr)
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)
library(lattice)
library(MASS)
library(lme4)
library(nlme)
library(emmeans)
library(gamm4)
library(dplyr)
library(magrittr)
library(reshape2)
library(mgcv)
library(lmerTest)
library(lubridate)
library(car)
library(lmerTest)
library(pls)

#loading workspace
setwd("/Users/lske0002/Projects/mycorrhizal removal/Response Curves/data")

### Growing Season Data for Point Measurements
# missing data in plots in 2017 were gap-filled with treatment averages for that day

data.gs <- read.csv("/Users/lske0002/Projects/mycorrhizal removal/Response Curves/data/Tsoil&VWC_fourdayavg_noNov_13JULY2020_treatmentavg_gapfill.csv") 

#growing season (no Nov data)
data.gs$Block <- as.factor(data.gs$Block) 
data.gs$erm_shrubs <- as.factor(data.gs$erm_shrubs)
data.gs$ecm_roots <- as.factor(data.gs$ecm_roots)
data.gs$Date <- as.Date(data.gs$Date)
data.gs$Month <- as.numeric(format(data.gs$Date, format = "%m"))
data.gs$d <- decimal_date(data.gs$Date) #making a vector of years with fractions
data.gs$Year<- as.factor(data.gs$Year)
data.gs$plot <- paste(data.gs$Treatment,data.gs$Block, sep = "")
data.gs$Respiration_mg <- (data.gs$Respiration)*(12.01/44.01)*1000

#subset by treatment
data.gs_DC <- data.gs[data.gs$Treatment == "DC",] #disturbance control, n = 176
data.gs_S <- data.gs[data.gs$Treatment == "S",] # combined removal, n = 176
data.gs_SE <- data.gs[data.gs$Treatment == "SE",] # pine root removal, n = 176
data.gs_SP <- data.gs[data.gs$Treatment == "SP",] # shrub removal, n = 176
data.gs_SPE <- data.gs[data.gs$Treatment == "SPE",] # control, n = 176

#combining treatments based on downstream analyses
data.gs_noDC <- rbind(data.gs_S,data.gs_SE,data.gs_SP,data.gs_SPE) # n = 704


######check to see if there are soil moisture and temperature differences from the 1-4 day averages before respiration measurements were taken

# lme treatment effects on Tsoil with Block as a random factor on growing season averages
# unbalanced design without filled in averages show very similar results

lme.tsoil <-lme(Tsoil_mean ~ erm_shrubs*ecm_roots+Year,
                random = ~ 1|Block/plot,
                data = data.gs_noDC, na.action = na.omit, correlation=corCAR1(form=~d|Block/plot))
anova(lme.tsoil)
plot(resid(lme.tsoil))
#write.csv(a,file = "lme.tsoil.point.msrmts.csv")

# lme treatment effects on VWC with Block as a random factor on growing season averages
# unbalanced data with less values in 2017 show a similar result although year is more significant

lme.vwc <-lme(VWC_mean ~ erm_shrubs*ecm_roots+Year, #year interaction with binary factors is not significant
              random = ~ 1|Block/plot,
              data = data.gs_noDC, na.action = na.omit, correlation=corCAR1(form=~d|Block/plot))
anova(lme.vwc)
coef(lme.vwc)
plot(resid(lme.vwc))
#write.csv(a,file = "lme.vwc.point.msrmts.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------#

# 1a full respiration model with interactions, no DC
lme1a.gs <- lme(log(Respiration) ~ ecm_roots*erm_shrubs*VWC_mean*Tsoil_mean + -(ecm_roots*erm_shrubs*VWC_mean*Tsoil_mean),
                random = ~ 1|Block/plot, 
                correlation=corCAR1(form=~d|Block/plot), # correlation structure of time series in nested design informed by decimal dates (d)
                data = na.omit(data.gs_noDC)) 

summary(lme1a.gs)

#checks for multicollinearity
vif(lme1a.gs)


lme1a.gs.anova <- anova(lme1a.gs)
lme1a.gs.coef <- coef(lme1a.gs)
lme1a.gs.anova

plot(lme1a.gs)
#write.csv(lme1a.coef,"/.csv")
#write.csv(lme1a.Anova,"") 

# 1b full respiration model with insignificant interactions and hypotheses, no DC
lme1b.gs <- lme(log(Respiration) ~ ecm_roots*VWC_mean*Tsoil_mean+erm_shrubs*VWC_mean*Tsoil_mean + ecm_roots*erm_shrubs,
                random = ~ 1|Block/plot, 
                correlation=corCAR1(form=~d|Block/plot), # correlation structure of time series in nested design informed by decimal dates (d)
                data = na.omit(data.gs_noDC)) 

summary(lme1b.gs)
anova(lme1b.gs)

#checks for multicollinearity
vif(lme1b.gs)

lme1b.gs.anova <- anova(lme1b.gs)
coef(lme1b.gs)

plot(lme1a.gs)
#write.csv(lme1a.coef,".csv")
#write.csv(lme1b.gs.anova,"lme1b.gs.anova.csv") 

##### Growing Season SPE vs. SPE-D

# growing season averaging soil respiration data across blocks by each date measurements were taken
data_grouped <- data.gs_SPEDC %>%
  group_by(Year,Date,Treatment) %>% 
  na.omit() %>%
  summarize_at(c("Respiration"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
data.gs_SPEDC_plot <- ungroup(data_grouped)


#colors
cbp2 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442")

#plot DC (SPE-D) and Control (SPE)
df <- data.frame(grp = data.gs_SPEDC_plot$Date,
                 fit = data.gs_SPEDC_plot$mean,
                 se = data.gs_SPEDC_plot$se)
j <- ggplot(df, aes(grp, fit, 
                    ymin = fit-se, 
                    ymax = fit+se, 
                    color=data.gs_SPEDC_plot$Treatment, 
                    shape = data.gs_SPEDC_plot$Treatment)) +
  geom_pointrange(size = .5) +
  geom_line(aes(linetype = data.gs_SPEDC_plot$Treatment))+
  scale_linetype_manual(values=c("dotted", "solid","solid"),
                        labels = c("Disturbance Control","Pine Root Removal","Control"),
                        name = "Treatment")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 18,angle = 60, vjust = .6),
        axis.text.y = element_text(size = 20),
        legend.justification = c("left","top"),
        legend.position = c(.05, .98),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  xlab(element_blank()) +
  ylab(expression(Respiration*" "*Rate*" "*mg*" "*C*" "*m^-2*" "*h^-1)) +
  scale_colour_manual(values = c("#D55E00","#009E73", "#E69F00" ),
                      labels = c("Disturbed Control","Pine Root Removal","Control"), 
                      name = "Treatment") +
  scale_shape_manual(values=c(7,19,19),
                     labels = c("Disturbed Control","Pine Root Removal","Control"), 
                     name = "Treatment") +
  scale_x_date(limits=start.end, 
               date_breaks = "2 months",
               date_labels = "%b")
j

year1.gs_SPEDC <- data.gs_SPEDC[data.gs_SPEDC$Year1 == "1",]
year2.gs_SPEDC <- data.gs_SPEDC[data.gs_SPEDC$Year1 == "2",]
year3.gs_SPEDC <- data.gs_SPEDC[data.gs_SPEDC$Year1 == "3",]

#DC stats
lme_dc.gs <- lme(log(Respiration) ~ Treatment,
                 random = ~ 1|Block/plot, 
                 correlation=corCAR1(form=~d|Block/plot), # correlation structure of time series in nested design informed by decimal dates (d)
                 data = na.omit(year3.gs_SPEDC))

anova(lme_dc.gs)

dc.emmeans <- emmeans(lme_dc.gs, c("Treatment"), type = "response")
pairs(dc.emmeans, type = "response", adjust = "bonf")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


######## Not Gap Filling Soil Moisture and Temperature from Missing Loggers  ######

#no gap filled averages for 2017
data.gs.miss <- read.csv("/Users/lske0002/Projects/mycorrhizal removal/Response Curves/data/Tsoil&VWC_fourdayavg_noNov_13JULY2020_missingLoggers.csv") 

str(data.gs.miss)

#growing season (no Nov data)
data.gs.miss$Block <- as.factor(data.gs.miss$Block) 
data.gs.miss$erm_shrubs <- as.factor(data.gs.miss$erm_shrubs)
data.gs.miss$ecm_roots <- as.factor(data.gs.miss$ecm_roots)
data.gs.miss$Date <- as.Date(data.gs.miss$Date)
data.gs.miss$Month <- as.numeric(format(data.gs.miss$Date, format = "%m"))
data.gs.miss$d <- decimal_date(data.gs.miss$Date) #making a vector of years with fractions
data.gs.miss$Year1<- ifelse(data.gs.miss$Year=="2017",1,ifelse(data.gs.miss$Year=="2018",2,3))
data.gs.miss$plot <- paste(data.gs.miss$Treatment,data.gs.miss$Block, sep = "")

#subset by treatment
data.gs.miss_DC <- data.gs.miss[data.gs.miss$Treatment == "DC",]
data.gs.miss_S <- data.gs.miss[data.gs.miss$Treatment == "S",] # n = 200
data.gs.miss_SE <- data.gs.miss[data.gs.miss$Treatment == "SE",] # n = 200
data.gs.miss_SP <- data.gs.miss[data.gs.miss$Treatment == "SP",] # n = 200
data.gs.miss_SPE <- data.gs.miss[data.gs.miss$Treatment == "SPE",] # n = 200

#combining treatments based on downstream analyses
#exclude the disturbance control (DC)
data.gs.miss_noDC <- rbind(data.gs_S,data.gs_SE,data.gs_SP,data.gs_SPE) 
#changing units from g CO2 to mg C
data.gs.miss_noDC$Respiration_mg <- (data.gs_noDC$Respiration)*(12.01/44.01)*1000

lme.vwc <-lme(VWC_mean ~ erm_shrubs*ecm_roots+Year, #year interaction with binary factors is not significant
              random = ~ 1|Block/plot,
              data = data.gs.miss_noDC, na.action = na.omit, correlation=corCAR1(form=~d|Block/plot))
anova(lme.vwc)
coef(lme.vwc)
plot(resid(lme.vwc))

data.gs.blocks257 <- subset(data.gs.miss_noDC, Block %in% c(2,5,7))


lme.vwc.257 <-lme(VWC_mean ~ erm_shrubs*ecm_roots+Year, #year interaction with binary factors is not significant
              random = ~ 1|Block/plot,
              data = data.gs.blocks257 , na.action = na.omit, correlation=corCAR1(form=~d|Block/plot))
anova(lme.vwc.257)
coef(lme.vwc.257)
plot(resid(lme.vwc.257))
