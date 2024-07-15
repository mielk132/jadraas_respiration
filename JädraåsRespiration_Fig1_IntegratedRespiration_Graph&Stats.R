# Jädraås Respiration Paper
#Figure 1 - Integrated Soil Respiration and Statistics

#May 25 2021
#Updated OCT 29 2021

rm(list=ls()) #clear the workspace


#libraries
library(mgcv); library(nlme); library(ggplot2);library(tidyverse)

setwd("/Users/lske0002/Projects/mycorrhizal removal/Respiration/data/Integrated")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#growing season

#mean values

growingseason.plot.co2 <- read.csv("Integrated Growing Season Plot means at Jädraås IhV reformatted.csv")
str(growingseason.plot.co2 <- read.csv("Integrated Growing Season Plot means at Jädraås IhV reformatted.csv"))
# averaging growing season integrated results for each plot and year
data_grouped <- growingseason.plot.co2 %>%
  group_by(Year,Treatment) %>% 
  na.omit() %>%
  summarize_at(c("g_C_m.2"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
integrated.gs_mean_by.year <- ungroup(data_grouped)
integrated.gs_mean_by.year <- as.data.frame(integrated.gs_mean_by.year)

#write.csv(integrated.gs_mean_by.year, "integrated.gs_mean_by.year.csv" )

# averaging growing season integrated results for each plot across all three years
data_grouped <- growingseason.plot.co2 %>%
  group_by(Treatment) %>% 
  na.omit() %>%
  summarize_at(c("g_C_m.2"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
integrated.gs_mean <- ungroup(data_grouped)
integrated.gs_mean
growingseason.plot.co2$Year <- as.factor(growingseason.plot.co2$Year)

#write.csv(integrated.gs_mean, "integrated.gs_mean.csv" )

#growing season (end of April or May-October) graph
# condition = year 1, 2, 3 

df <- data.frame(grp = integrated.gs_mean_by.year$Treatment, fit = integrated.gs_mean_by.year$mean, se = integrated.gs_mean_by.year$se, condition = as.factor(integrated.gs_mean_by.year$Year))

l <- ggplot(df, aes(fill=condition, y=fit, x=grp)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  theme_classic() +
  geom_errorbar(aes(ymax=fit+se, ymin = fit-se),
                position=position_dodge(.9),width=.2)+
  ggtitle("Jädraås IhV Respiration Growing Season 2017-2019") +
  xlab(element_blank()) +
  ylab(expression("g C m-2 yr-1"))+
  scale_y_continuous(expand = (c(0,0)))+
  scale_x_discrete(limits = c("Control", "Shrub Removal", "Pine Root Exclusion","Shrub & Pine Root Exclusion"))
l

#stats

#log() to satisfy normal distribution of residuals - it does not change the interpretation. 

lme1.gs <- lme(log(g_C_m.2) ~ EcM_root*ErM_root*Year, 
               random = ~ 1|plot/Block,
               data = growingseason.plot.co2)
summary(lme1.gs)
anova(lme1.gs)
plot(resid(lme1.gs))

#The interaction with year and ecm roots was taken out (P > 0.2) because there were no significant effects but kept in ecm*erm interaction to test the hypothesis whether there is some sort of competitive release. 

lme2.gs <- lme(log(g_C_m.2) ~ EcM_root*ErM_root+ErM_root*Year+EcM_root*Year, 
               random = ~ 1|plot/Block,
               data = growingseason.plot.co2)
summary(lme2.gs)
integrated.anova <- anova(lme2.gs)
plot(resid(lme2.gs))
coef(lme2.gs)
integrated.anova

write.csv(integrated.anova, file = "integrated.anova.csv")
