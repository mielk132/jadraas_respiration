##Nitrogen Availability from Resin Strips 2017-2018
##Lingon Leaf Graphs IRMS 2017

#Author: Louis Mielke
#Date Created: 2020 March 3 
#Date Updated: 2021 May 31

#testing the hypothesis that lingon leaves would increase access to nitrogen sources with the decay of roots and mycelium

rm(list=ls()) #clear the workspace

#libraries
library(mgcv); library(nlme); library(ggplot2); library(multcomp); library(car); library(ape)
library(itsadug);library(MASS);library(grid);library(readxl);library(ggplot2); library(readxl);
library(ggsci);library(tidyr);library(tidyverse);library(scales);library(gridExtra)
library(grid);library(lattice);library(MASS);library(lme4);library(nlme);library(emmeans);library(gamm4);library(dplyr);library(magrittr);library(reshape2)
library(mgcv);library(lmerTest);library(lubridate);library(car);library(lmerTest);library(pls); library(forcats)

setwd("/Users/lske0002/Projects/mycorrhizal removal/Resin Strips")
n_data <- read.csv("ResinStrips_2years.csv")
n_data$NH4.N_mg.L <- n_data$NH4.N_ug.L/1000
n_data$NH4.N_mg.L[n_data$NH4.N_mg.L < .02] = 0

#summary of n data by year, season and Treatment from Block seasonal averages (n=8)
n_data_avg <- na.omit(n_data) %>%
  group_by(box, Treatment) %>%
  drop_na() %>% 
  summarize_at(c("NH4.N_mg.L","NO3.N_ug.L"), funs(mean, sd, n()))

n_data_avg
#write.csv(n_data_avg,file = "jädraås_n_data_avg.csv")

#nitrate levels are so low that we are not going to continue analyzing them

#setting start and end date for each sample as a recognized date
start <- as.Date(n_data$start_date, format = "%Y-%m-%d")
end <- as.Date(n_data$end_date, format = "%Y-%m-%d")
days <- difftime(end, start, units = "days") #the number of days incubated
n_data$days <- days

#In 2017, 49 day incubation 19-05-2017 to 07-07-2017 (period 1) and 63 day incubation 07-07-2017 to 08-09-2017 (period 2)
#In 2018, 118 day incubation 22-05-2018 to 17-09-2018 (period 1) and 59 day incubation 17-09-2018 to 15-11-2018 (period 2)

#subsetting datasets into years
n_data2017 <- n_data[n_data$Year == "2017",]
n_data2018 <- n_data[n_data$Year == "2018",]
n_data2017_spring <- subset(n_data2017, season %in% c("spring"))
n_data2017_summer <- subset(n_data2017, season %in% c("summer"))
n_data2018_fall <- subset(n_data2018, season %in% c("fall"))
n_data2018_summer <- subset(n_data2018, season %in% c("summer"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

par(mfrow = c(1,1)) #setting the graphing window in base R


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#2017 boxplot of NH4 concentrations for all treatments
NH4.box2017 <- ggplot(data = na.omit(n_data2017)) +
  aes(x = Treatment, y = NH4.N_mg.L, fill = Treatment) +
  geom_boxplot() +
  theme_bw()+
  ylim(0,50)+
  facet_wrap(vars(season))
NH4.box2017

#linear mixed model (two-way ANOVA) testing treatment effects for 2017
#sqrt satisfies variation in residuals
lmeNH4.2017 <- lme(sqrt(NH4.N_mg.L) ~ Treatment*season, random = ~ 1|Block, data = na.omit(n_data2017)) 
anova(lmeNH4.2017)
plot(resid(lmeNH4.2017))

emmeans.NH417 <- emmeans(lmeNH4.2017, c("Treatment","season"), type = "response")
pairs(emmeans.NH417, type = "response", adjust = "bonf")
contrast(regrid(emmeans.NH417))

#one-way ANOVAs for sig letters for the bar graph

#spring 2017
lmeNH4.2017.spring <- lme(sqrt(NH4.N_mg.L) ~ Treatment, random = ~ 1|Block, data = na.omit(n_data2017_spring)) 
anova(lmeNH4.2017.spring)
plot(resid(lmeNH4.2017.spring))

emmeans.NH417.spring <- emmeans(lmeNH4.2017.spring, c("Treatment"), type = "response")
pairs(emmeans.NH417.spring, type = "response", adjust = "bonf")
contrast(regrid(emmeans.NH417))

#summer 2017
lmeNH4.2017.summer <- lme(sqrt(NH4.N_mg.L) ~ Treatment, random = ~ 1|Block, data = na.omit(n_data2017_summer)) 
anova(lmeNH4.2017.summer)
plot(resid(lmeNH4.2017.summer))

emmeans.NH417.summer <- emmeans(lmeNH4.2017.summer, c("Treatment"), type = "response")
pairs(emmeans.NH417.summer, type = "response", adjust = "bonf")
contrast(regrid(emmeans.NH417.summer))

#linear mixed model testing treatment effects for 2017
#sqrt satisfies variation in residuals
lmeNH4.2017.factorial <- lme(sqrt(NH4.N_mg.L) ~ shrub*pine+pine*season+shrub*season, random = ~ 1|Block, data = na.omit(n_data2017)) 
anova(lmeNH4.2017.factorial)
coef(lmeNH4.2017.factorial)
plot(resid(lmeNH4.2017.factorial))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#2018 boxplot of NH4 concentrations for all treatments
NH4.box2018 <- ggplot(data = na.omit(n_data2018)) +
  aes(x = Treatment, y = NH4.N_mg.L, fill = Treatment) +
  geom_boxplot() +
  theme_bw()+
  facet_wrap(vars(season))
NH4.box2018

#linear mixed model testing treatment effects for 2018
#sqrt satisfies variation in residuals
lmeNH4.2018 <- lme(sqrt(NH4.N_mg.L) ~ Treatment*season, random = ~ 1|Block, data = na.omit(n_data2018))
anova(lmeNH4.2018)
plot(resid(lmeNH4.2018))

emmeans.NH418 <- emmeans(lmeNH4.2018, c("Treatment","season"), type = "response")
pairs(emmeans.NH418, type = "response", adjust = "bonf")

#one-way ANOVAs for the bar graph

#summer season
lmeNH4.2018.summer <- lme(log(NH4.N_mg.L) ~ Treatment, random = ~ 1|Block, data = na.omit(n_data2018_summer))
anova(lmeNH4.2018.summer)
plot(resid(lmeNH4.2018.summer)
     
#post-hoc test
emmeans.NH418.summer <- emmeans(lmeNH4.2018.summer, c("Treatment"), type = "response")
pairs(emmeans.NH418.summer, type = "response", adjust = "bonf")

#fall season
lmeNH4.2018.fall <- lme(sqrt(NH4.N_mg.L) ~ Treatment, random = ~ 1|Block, data = na.omit(n_data2018_fall))
anova(lmeNH4.2018.fall)
plot(resid(lmeNH4.2018.fall))

#post-hoc test
emmeans.NH418.fall <- emmeans(lmeNH4.2018.fall, c("Treatment"), type = "response")
pairs(emmeans.NH418.fall, type = "response", adjust = "bonf")


#factorial linear mixed model testing treatment effects for 2018
#sqrt satisfies variation in residuals
lmeNH4.2018 <- lme(sqrt(NH4.N_mg.L) ~ pine*shrub+pine*season, random = ~ 1|Block, data = na.omit(n_data2018))
anova(lmeNH4.2018)
plot(resid(lmeNH4.2018))
#------------------------------
n_data_avg$Treatment <- factor(n_data_avg$Treatment, levels = c("Control","Disturbed Control","Shrub Removal","Pine Root Removal","Combined Removal"))


n_data_avg$NH4_se <- (n_data_avg$NH4.N_mg.L_sd)/sqrt(n_data_avg$NH4.N_mg.L_n) #calculating standard error

cbp2 <- c("#E69F00","#D55E00", "#56B4E9", "#009E73",
          "#F0E442")

df <- data.frame(Treatment = n_data_avg$Treatment, fit = n_data_avg$NH4.N_mg.L_mean, condition = n_data_avg$box)

nh4 <- ggplot(df, aes(fill=Treatment, y=fit, x=condition)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20),
        legend.justification = c("left","top"),
        legend.position = c(.05, .85),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  geom_errorbar(aes(ymax=fit+n_data_avg$NH4_se, ymin = fit-n_data_avg$NH4_se),
                position=position_dodge(.9),width=.2)+
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"))+
  xlab("Year 1                         Year 2") +
  ylab(expression(mg*" "*NH["4"]^"+"*"- N"*"  "*L^-1))+
  scale_y_continuous(limits = c(0,45),expand = (c(0,0)))+
  scale_x_discrete(limits = c("a","b", "c", "d"), labels = c("Spring 17","Summer 17","Summer 18","Fall 18"))+
  annotate("text", x=.75, y=43, label= "(a)", size = 10)
nh4

### Lingon Leaf Isotopes

#loading lingon leaf IRMS data
setwd("/Users/lske0002/Projects/mycorrhizal removal/LingonLeaves/data")
lingon.raw <- read.csv("Lingon_Leaves_2017_clean.csv")


lingon.raw$Treatment <- as.factor(as.character(lingon.raw$Treatment))

# averaging lingon leave data across blocks by each measurement

data_grouped <- lingon.raw %>%
group_by(Treatment) %>%
drop_na() %>% 
summarize_at(c("N","delta15N","C","delta13C","CN"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
lingon <- ungroup(data_grouped)


#plot %C
df <- data.frame(grp = lingon$Treatment, fit = lingon$C_mean, se = lingon$C_se)
a <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=lingon$Treatment)) +
geom_pointrange(size = 1) +
theme_classic() +
ggtitle("     Lingon leaves %C across Treatments") +
xlab(element_blank()) +
ylab("%C") +
scale_colour_manual(values = c("black","black","black","black","black")) +
theme(text = element_text(size = 15) , 
axis.line = element_line( size = 1, linetype = "solid"), 
axis.text=element_text(size=19),
legend.position = "none")
a

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#dotplot %N
df <- data.frame(grp = lingon$Treatment, fit = lingon$N_mean, se = lingon$N_se)
b <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=lingon$Treatment)) +
geom_pointrange(size = 1) +
theme_classic() +
ggtitle("         Lingon leaves %N in sept 2017") +
xlab(element_blank()) +
ylab("%N") +
scale_colour_manual(values = c("black","black","black","black","black")) +
theme(text = element_text(size = 15) , 
axis.line = element_line( size = 1, linetype = "solid"), 
axis.text=element_text(size=19),
legend.position = "none")
b

#barplot %N

lingon$Treatment <- factor(lingon$Treatment, levels = c("Control","Disturbed Control","Shrub Removal","Pine Root Removal","Shrub & Pine Root Removal"))

df <- data.frame(grp = lingon$Treatment, fit = lingon$N_mean, se = lingon$N_se)
leaf_n <- ggplot(df, aes(fill=grp, y=fit, x=grp)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20),
        legend.position = "none",
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  geom_errorbar(aes(ymax=fit+se, ymin = fit-se),
                position=position_dodge(.9),width=.2)+
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"))+
  ylab(expression("leaf N concentration (% of dry mass)"))+
  scale_y_continuous(limits = c(0,2.1),expand = (c(0,0)))+
  annotate("text", x=.75, y=2, label= "(b)", size = 10)
leaf_n


#stats %N
lme_N <- lme(N ~ Treatment, random = ~ 1|Block, data = na.omit(lingon.raw))
anova(lme_N)
plot(resid(lme_N))
emmeans.N <- emmeans(lme_N, c("Treatment"), type = "response")
pairs(emmeans.N, type = "response", adjust = "bonf")


#stats %N
lme_N <- lme(N ~ pine*shrub, random = ~ 1|Block, data = na.omit(lingon.raw))
anova(lme_N)
plot(resid(lme_N))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#plot C:N
df <- data.frame(grp = lingon$Treatment, fit = lingon$CN_mean, se = lingon$CN_se)
c <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=lingon$Treatment)) +
geom_pointrange(size = 1) +
theme_classic() +
ggtitle("     Lingon leaves C:N across Treatments") +
xlab(element_blank()) +
ylab("C:N") +
scale_colour_manual(values = c("black","black","black","black","black")) +
theme(text = element_text(size = 15) , 
axis.line = element_line( size = 1, linetype = "solid"), 
axis.text=element_text(size=19),
legend.position = "none")
c

#stats C:N
lme_CN <- lme(CN ~ Treatment, random = ~ 1|Block, data = na.omit(lingon.raw)) 
anova(lme_CN)
plot(resid(lme_CN))
emmeans.CN <- emmeans(lme_CN, c("Treatment"), type = "response")
pairs(emmeans.CN, type = "response", adjust = "bonf")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#plot delta15N

plot(lingon.raw$Block,lingon.raw$delta15N) 

df <- data.frame(grp = lingon$Treatment, fit = lingon$delta15N_mean, se = lingon$delta15N_se)
d <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=lingon$Treatment)) +
geom_pointrange(size = 1) +
theme_classic() +
ggtitle("     Lingon leaves δ¹⁵N (‰) across Treatments") +
xlab(element_blank()) +
ylab("δ¹⁵N (‰)") +
scale_colour_manual(values = c("black","black","black","black","black")) +
theme(text = element_text(size = 15) , 
axis.line = element_line( size = 1, linetype = "solid"), 
axis.text=element_text(size=19),
legend.position = "none")
d
# delta15N stats 
lme_15N <- lme(delta15N ~ Treatment, random = ~ 1|Block, data = na.omit(lingon.raw)) 
anova(lme_15N)
plot(resid(lme_15N))
emmeans.15N <- emmeans(lme_15N, c("Treatment"), type = "response")
pairs(emmeans.15N, type = "response", adjust = "bonf")

# delta15N and %N by treatment

df <- data.frame(grp = lingon$N_mean, fit = lingon$delta15N_mean)

n15n <- ggplot(df, 
aes(grp, fit, color=lingon$Treatment, shape = lingon$Treatment))+
geom_point(aes(size = 14))+
geom_errorbar(aes(ymin = fit+lingon$delta15N_se,ymax = fit-lingon$delta15N_se))+
geom_errorbarh(aes(xmin = grp+lingon$N_se,xmax = grp-lingon$N_se))+

theme(panel.border= element_rect(colour = "black", size=1, fill=NA),
panel.background = element_rect(colour = "black", size=1, fill=NA)) +
xlab("%N") +
ylab(expression(paste(delta^{15}, "N \u2030"))) +
scale_shape_manual(values=c(19,19,19,19,19)) +
scale_color_manual(values = c(c))+
theme(legend.position = "none",
axis.line = element_line( size = 1, linetype = "solid"), 
axis.title=element_text(size=16,face="bold"),
axis.text.x = element_text(size = 16),
axis.text.y = element_text(size = 20))+
annotate("text", x=1.05, y=-3.5, label= "(b)", size = 10)
n15n


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#plot delta13C
df <- data.frame(grp = lingon$Treatment, fit = lingon$delta13C_mean, se = lingon$delta13C_se)
e <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=lingon$Treatment)) +
geom_pointrange(size = 1) +
theme_classic() +
ggtitle("     Lingon leaves δ¹³C (‰) across Treatments") +
xlab(element_blank()) +
ylab("δ¹³C (‰)") +
scale_colour_manual(values = c("black","black","black","black","black")) +
theme(text = element_text(size = 15) , 
axis.line = element_line( size = 1, linetype = "solid"), 
axis.text=element_text(size=19),
legend.position = "none")
e

#stats delta13C
lme_13C <- lme(delta13C ~ Treatment, random = ~ 1|Block, data = na.omit(lingon.raw)) 
anova(lme_13C)
plot(resid(lme_13C))
emmeans.13C <- emmeans(lme_13C, c("Treatment"), type = "response")
pairs(emmeans.13C, type = "response", adjust = "bonf")


## combining both NH4 and %N graphs

require(gridExtra)
plot1 <- nh4
plot2 <- leaf_n
grid.arrange(plot1, plot2, ncol=2)

setwd("/Users/lske0002/Projects/mycorrhizal removal/Writing/Respiration Paper/final graphs and tables")

# Open a pdf file
pdf("rplot_NH4_leafN.pdf") 

grid.arrange(plot1, plot2, ncol=2)

# Close the pdf file
dev.off()
