# "Jädraås Respiration Paper Figure 2a-c
# 
# Supplementary Table 1 & Supp Fig 4
#Created MAY 3 2021
#Updated 9 NOV 2021

rm(list=ls()) #clear workspace


#libraries
library(lubridate);library(tidyverse);library(ggplot2);library(dplyr
library(mgcv); library(nlme); library(ggplot2);library(car); library(ape);library(readxl);library(decimaldate);library(MASS);library(grid);library(readxl);library(ggsci);library(tidyr);;library(scales);library(gridExtra);library(grid);library(lattice);library(MASS);library(lme4);library(nlme);library(emmeans);library(gamm4);library(magrittr);library(reshape2);library(mgcv);library(lmerTest);library(car);library(lmerTest);library(pls); library(forcats);library(patchwork)

####### Whole Exprmnt Averages of Moisture and Temp #####
# loading all soil temp and vwc data

#setwd("/")

#vwc & temp loggers 2016Dec07-2019Nov20 long format
env.data <- read.csv("Data_1.csv") 
env.data <- na.omit(env.data)

#characterizing time categories
env.data$Date <- as.Date(env.data$Date)
env.data$Treatment <- factor(env.data$Treatment, levels = c('C', 'E', 'T','TE','DC'))
#add month as categorical variable
env.data$month <- as.factor(as.numeric((format(env.data$Date, format = "%m"))))

env.data$d <- as.numeric(decimal_date(env.data$Date)) # making a vector of years with fractions

# setting negative values of moisture to 0
env.data$VWC <- ifelse(env.data$VWC < 0, 0, env.data$VWC) 

#making Block a factor, not an integer & setting up a column with plot ID
env.data$Block <- as.factor(env.data$Block)
env.data$plot <- as.factor(paste(env.data$Treatment,env.data$Block))  

#summary for vwc data for each day over all plots
vwc_jädraås1 <- env.data %>%
  group_by(Date) %>% 
  drop_na() %>% 
  summarize_at(c("VWC"), funs(mean,n(),se=sd(.)/sqrt(n())))
vwc_jädraås1

vwc_jädraås1$Treatment <- rep("All",1079) #dummy variables


# Define Start and end times for the subset as R objects that are the time class
startTime <- as.Date("2017-01-01")
endTime <- as.Date("2019-12-01")

start.end <- c(startTime,endTime)
start.end

#plot soil moisture full data series average (daily)
df <- data.frame(date = vwc_jädraås1$Date, fit_vwc_avg = vwc_jädraås1$mean)
vwc_daily_average <- ggplot(df, aes(date, fit_vwc_avg)) +
  geom_line()+
  theme_classic() +
  ggtitle("     ") +
  xlab(element_blank())+
  ylab(expression("Soil moisture")) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b")+
  scale_y_continuous(limits = c(0,.35),expand = c(0.00,0.00))+
  theme(text = element_text(size = 12) , 
        axis.line = element_line( size = .5, linetype = "solid"),axis.text.y = element_text(size = 16),
        panel.border= element_rect(colour = "black", size=.5, fill=NA),
        panel.background = element_rect(colour = "black", size=.5, fill=NA))+
  scale_color_manual(values = c("#737373"))
vwc_daily_average

#summary for tsoil data for each day over all plots
tsoil_jädraås1 <- env.data %>%
  group_by(Date) %>% 
  drop_na() %>% 
  summarize_at(c("Tsoil"), funs(mean,n(),se=sd(.)/sqrt(n())))
tsoil_jädraås1

tsoil_jädraås1$Treatment <- rep("All",1079) #add dummy variables
tsoil_jädraås1$Year <- rep("NA",1079) #add dummy variables

#plot soil temp full data series average (daily) 
df <- data.frame(grp = tsoil_jädraås1$Date, fit = tsoil_jädraås1$mean)
tsoil_daily_average <- ggplot(df, aes(grp, fit)) +
  geom_line()+
  theme_classic() +
  ggtitle("     ") +
  xlab(element_blank())+
  ylab(expression("Soil temperature")) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b")+
  scale_y_continuous(limits = c(0,22),expand = c(0.00,0.00))+
  theme(text = element_text(size = 12) , 
        axis.line = element_line( size = .5, linetype = "solid"),axis.text.y = element_text(size = 16),
        panel.border= element_rect(colour = "black", size=.5, fill=NA),
        panel.background = element_rect(colour = "black", size=.5, fill=NA))+
  scale_color_manual(values = c("#737373"))
tsoil_daily_average

######### Figure 2a Respiration Point Measurements  ######
data.gs <- read.csv("Data_2.csv") 

#growing season (no November data)
data.gs$Block <- as.factor(data.gs$Block) 
data.gs$erm_shrubs <- as.factor(data.gs$shrubs) 
data.gs$ecm_roots<- as.factor(data.gs$pine_roots)
data.gs$Date <- as.Date(data.gs$Date)
data.gs$Month <- as.numeric(format(data.gs$Date, format = "%m"))
data.gs$d <- decimal_date(data.gs$Date) #making a vector of years with fractions
data.gs$Year1<- ifelse(data.gs$Year=="2017",1,ifelse(data.gs$Year=="2018",2,3))
data.gs$plot <- paste(data.gs$Treatment,data.gs$Block, sep = "")

#subset by treatment
data.gs_DC <- data.gs[data.gs$Treatment == "DC",]
data.gs_TE <- data.gs[data.gs$Treatment == "TE",] # n = 200
data.gs_T <- data.gs[data.gs$Treatment == "T",] # n = 200
data.gs_E <- data.gs[data.gs$Treatment == "E",] # n = 200
data.gs_C <- data.gs[data.gs$Treatment == "C",] # n = 200

#combining treatments based on downstream analyses
data.gs_noDC <- rbind(data.gs_TE,data.gs_T,data.gs_E,data.gs_C)
#changing units from (g CO2 m^2) to (mg C m^2)
data.gs$Respiration_mg <- (data.gs$Respiration)*(12.01/44.01)*1000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# averaging soil respiration data across blocks by each date measurements were taken
data_grouped <- resp.gs %>%
  group_by(Year,Date,Treatment) %>%
  drop_na() %>%
  summarize_at(c("Respiration_mg"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
resp.gs.avg <- ungroup(data_grouped)

# Define Start and end times for the subset as R objects that are the time class
startTime <- as.Date("2017-01-1")
endTime <- as.Date("2019-12-01")

start.end <- c(startTime,endTime)
start.end

#colors for DC, 
cbp2 <- c("#E69F00", "black","#56B4E9","#009E73","#E6CF02")
lines <- c("solid","dashed","solid","solid","solid")

#plot growing season Jädraås mycorrhizal roots and saprotrophic respiration rates 2017-2019
df <- data.frame(grp = resp.gs.avg$Date, fit = resp.gs.avg$mean, se = resp.gs.avg$se, treatment = resp.gs.avg$Treatment)
figure2a_resp <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color= treatment, linetype = treatment)) +
  xlab(element_blank()) +
  ylab(expression(Soil*" "*respiration*" "*mg*" "*C*" "*m^-2*" "*h^-1)) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b") +
  scale_y_continuous(expand =c(0,0),limits = c(0,415)) +
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        #legend.title = element_text(size = 12),
        #legend.text = element_text(size = 12),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=.5, fill=NA)) +
  geom_line() +
  #geom_errorbar(data = df[df$treatment == "DC",], linetype = "solid")+
  scale_linetype_manual(values = c(lines))+
  geom_pointrange(data = df[df$treatment != "DC",], size = 1, shape =1) +
  scale_color_manual(values = c(cbp2)) +
  guides(color=guide_legend(title="Treatment"))
figure2a_resp


############## FIGURE 2b Soil Moisture Point Measurements ##########

# averaging soil moisture by date and treatment
data.vwc.point.avg <- data.gs_noDC %>%
  group_by(Date, Treatment) %>%
  drop_na() %>% 
  summarize_at(c("VWC_mean"), funs(mean, n(), se=sd(.)/sqrt(n())))
data.vwc.point.avg <- ungroup(data.vwc.point.avg)

#adding the all plot average of soil moisture (a line) to the data set with the point measurements of soil moisture

alla_moist <- rbind(data.vwc.point.avg,vwc_jädraås1)

#colors
cbp <- c("black","#E6CF02", "#009E73","#56B4E9",
         "#E69F00")
#plot 
df1 <- data.frame(grp = alla_moist$Date, fit = alla_moist$mean, se = alla_moist$se, treatment = alla_moist$Treatment)
figure2b_moist <- ggplot(df1, aes(shape = "1", grp, fit, ymin = fit-se, ymax = fit+se, color=treatment)) +
  theme_classic() +
  ggtitle("   ") +
  xlab(element_blank()) +
  ylab(expression("Soil moisture")) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b")+
  scale_y_continuous(limits = c(0,.35),expand = c(0.00,0.00))+
  theme(text = element_text(size = 12) , 
        axis.line = element_line( size = .5, linetype = "solid"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        panel.border= element_rect(colour = "black", size=.5, fill=NA),
        legend.position = "none",
        panel.background = element_rect(colour = "black", size=.5, fill=NA))+
  geom_line(data = df1[df1$treatment == "All",])+
  geom_pointrange(data = df1[df1$treatment != "All",],size = 1) +
  scale_color_manual(values = c(cbp))+
  scale_shape_manual(values=c(1,1,1,1))+
scale_fill_discrete(guide = guide_legend(reverse=TRUE))
figure2b_moist

# averaging soil moisture by year
data.vwc.point.avg <- data.gs_noDC %>%
  group_by(Year) %>%
  drop_na() %>% 
  summarize_at(c("VWC_mean"), funs(mean, n(), se=sd(.)/sqrt(n())))
data.vwc.point.avg <- ungroup(data.vwc.point.avg)


############# FIGURE 2c Soil Temp Point Measurements ###############

# averaging soil temp across blocks within treatments, dates and years
data.tsoil.point.avg <- data.gs_noDC %>%
  group_by(Treatment, Year, Date) %>%
  drop_na() %>% 
  summarize_at(c("Tsoil_mean"), funs(mean, n(), se=sd(.)/sqrt(n())))
data.tsoil.point.avg <- ungroup(data.tsoil.point.avg)


#adding the all plot average of soil temp to fill out the plot
alla_tsoil <- rbind(data.tsoil.point.avg,tsoil_jädraås1)

#plot 
df <- data.frame(grp = alla_tsoil$Date, fit = alla_tsoil$mean, se = alla_tsoil$se, treatment = alla_tsoil$Treatment)
figure2c_temp <- ggplot(df, aes(shape = "1", grp, fit, ymin = fit-se, ymax = fit+se, color=treatment)) +
  theme_classic() +
  ggtitle("     ") +
  xlab(element_blank()) +
  ylab(expression("Soil temperature")) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b")+
  scale_y_continuous(limits = c(-5,22),expand = c(0.00,0.00))+
  theme(text = element_text(size = 12) , 
        axis.line = element_line( size = .5, linetype = "solid"),
        axis.text.y = element_text(size = 16),
        panel.border= element_rect(colour = "black", size=.5, fill=NA),
        legend.position = "none",
        panel.background = element_rect(colour = "black", size=.5, fill=NA))+
  geom_line(data = df[df$treatment == "All",])+
  geom_pointrange(data = df[df$treatment != "All",],size = 1) +
  scale_color_manual(values = c(cbp))+
  scale_shape_manual(values=c(1,1,1,1))+
  scale_fill_discrete(guide = guide_legend(reverse=TRUE))
figure2c_temp

#combining plots into display
require(gridExtra)
grid.arrange(figure2a_resp, figure2b_moist, figure2c_temp, nrow=3, padding = c(0,0))

     
######## SUPPLEMENTARY FIGURE 4 ######

### ORGANIZING ###

#subset the dataset by growing season months (may-october)
env.data.gs <- subset(env.data, month %in% c("5","6","7","8","9","10"))

#add year as a categorical variable
env.data.gs$year <- as.factor(as.factor((format(env.data.gs$Date, format = "%y"))))

#subset the dataset by specific blocks (may-october)
env.data.gs.blocks257 <- subset(env.data.gs, Block %in% c(2,5,7))
env.data.gs.block5 <- subset(env.data.gs, Block %in% c(5))


#summary of jädraås climate data by year and Treatment from daily averages
env_data_gs_byplot <- env.data.gs %>%
  group_by(Treatment,Block,year) %>%
  drop_na() %>% 
  summarize_at(c("VWC","Tsoil"), funs(mean, min, max,sd, n()))
env_data_gs_byplot

#summary of jädraås climate data by year and Treatment from blocks 2,5,7 daily averages
env_data_gs_byplot_blocks257 <- env.data.gs.blocks257 %>%
  group_by(year, Treatment, Block) %>%
  summarize_at(c("VWC","Tsoil"), funs(mean, min, max,sd, n()))
env_data_gs_byplot_blocks257

#summary of jädraås climate data by year and Treatment from block 5 daily averages
#this is the only block in 2017 that had consistent soil temp and vwc for the combined shrub removal and pine root exclusion

env_data_gs_byplot_block5 <- env.data.gs.block5 %>%
  group_by(year, Treatment, Block) %>%
  summarize_at(c("VWC","Tsoil"), funs(mean, min, max,sd, n()))
env_data_gs_byplot_block5

env_data_gs_byplot_block5_TE <- env_data_gs_byplot_block5[env_data_gs_byplot_block5$Treatment == "TE",]
env_data_gs_byplot_blocks257_noTE <- env_data_gs_byplot_blocks257[env_data_gs_byplot_blocks257$Treatment %in% c("C","E","T","DC"),]

env_data_gs_17 <- subset(rbind(env_data_gs_byplot_block5_TE,env_data_gs_byplot_blocks257_noTE),year %in% c("17"))
env_data_gs_18_19 <- subset(env_data_gs_byplot, year %in% c("18","19"))


new_gs_data <- rbind(env_data_gs_17,env_data_gs_18_19)

write.csv(new_gs_data, file = "new_gs_data_2021NOV9.csv")

#new_gs_data <- read.csv("env_data_growingseason_new_noDC_factorial.csv")
#taking out averages with many missing dates

new_gs_data <- new_gs_data[new_gs_data$Tsoil_n > 3200,] 
str(new_gs_data)

#making averages across blocks for each treatment within each year
new_gs_data_avg <- new_gs_data %>%
  group_by(year, Treatment) %>%
  summarize_at(c("VWC_mean","Tsoil_mean"), funs(mean, min, max,sd, n()))
new_gs_data_avg

new_gs_data_avg$year <- as.numeric(new_gs_data_avg$year)
#new_gs_data_avg$year <- as.factor(new_gs_data_avg$year)
new_gs_data_avg$Tsoil_mean_se <- new_gs_data_avg$Tsoil_mean_sd / sqrt(new_gs_data_avg$Tsoil_mean_n)

#standard error
new_gs_data_avg$VWC_mean_se <- new_gs_data_avg$VWC_mean_sd / sqrt(new_gs_data_avg$VWC_mean_n)  

### GRAPHING ###

#plot soil temperature by growing season per year (n=8)
temp.jädraås.growingseason.avg <- ggplot(data = new_gs_data_avg) +
  aes(x = year, y = Tsoil_mean_mean, fill = Treatment) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Tsoil_mean_mean-Tsoil_mean_se,
                    ymax=Tsoil_mean_mean+Tsoil_mean_se), 
                width=.2, position=position_dodge(.9))+
  scale_y_continuous(limits = c(0,17),expand = (c(0,0)))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.justification = c("left","top"),
        legend.position = c("none"),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1, fill=NA))+
  ylab("soil temperature (°C)")+
  xlab("")+
  annotate("text", x=3.25, y=16.5, label= "(a)", size = 8)+
  scale_x_discrete(limits = c(1,2,3), labels = c("2017","2018","2019"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00"),labels = c("Control","Shrub Removal","Pine Root Removal","Pine Root & Shrub Removal", "Disturbed Control"))
temp.jädraås.growingseason.avg

#plot soil moisture by growing season per year (n=8)
vwc.jädraås.growingseason.avg <- ggplot(data = new_gs_data_avg) +
  aes(x = year, y = VWC_mean_mean, fill = Treatment) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=VWC_mean_mean-VWC_mean_se,
                    ymax=VWC_mean_mean+VWC_mean_se), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(limits = c(0,.25),expand = (c(0,0)))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.justification = c("right"),
        legend.position = c(.95, .85),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=1,fill=NA))+
  ylab(expression(paste("soil moisture (",m^3,"/",m^3,")",sep="")))+
  xlab("")+
  annotate("text", x=3.25, y=.26, label= "(b)", size = 8)+
  scale_x_discrete(limits = c(1,2,3), labels = c("2017","2018","2019"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00"),labels = c("Control","Shrub Removal","Pine Root Exclusion","Pine Root Exclusion & Shrub Removal", "Disturbed Control"))
vwc.jädraås.growingseason.avg 

require(gridExtra)
plot1 <- temp.jädraås.growingseason.avg
plot2 <- vwc.jädraås.growingseason.avg
grid.arrange(plot1, plot2, ncol=2)

#### STATISTICS for SUPP FIG 4 ######

str(env_data_gs_18_19) # only looking at year 2 and 3 because of missing data in the first year

env_data_gs_18_19$Block <- as.numeric(as.factor(env_data_gs_18_19$Block))
env_data_gs_18_19$plot <- as.factor(paste(env_data_gs_18_19$Treatment,env_data_gs_18_19$Block, sep="", collapse=NULL))
env_data_gs_18_19 <- env_data_gs_18_19[env_data_gs_18_19$Treatment != "DC",]
env_data_gs_18_19$year <- as.factor(env_data_gs_18_19$year)

# selecting averages with a certain criterion of hourly time points to make sure it is(3900)
env_data_gs_18_19 <- env_data_gs_18_19[env_data_gs_18_19$Tsoil_n > 3900,] 
str(env_data_gs_18_19)

library(nlme)

# lme treatment effects on Tsoil with Block as a random factor on growing season averages per treatment and year

lme.tsoil.18 <-lme(Tsoil_mean ~ Treatment,
                random = ~ 1|plot/Block,
                data = na.omit(subset(env_data_gs_18_19, year %in% c("18"))))
anova.t18 <- anova(lme.tsoil.18)
plot(resid(lme.tsoil.18))
qqnorm(resid(lme.tsoil.18))

library(emmeans)

#post hoc test
emmeans.tsoil.18 <- emmeans(lme.tsoil.18, c("Treatment"), type = "response")
pairs(emmeans.tsoil.18, type = "response", adjust = "bonf")


lme.tsoil.19 <-lme(Tsoil_mean ~ Treatment,
                   random = ~ 1|plot/Block,
                   data = na.omit(subset(env_data_gs_18_19, year %in% c("19"))))
anova.t19 <- anova(lme.tsoil.19)
plot(resid(lme.tsoil.19))
qqnorm(resid(lme.tsoil.19))

#post hoc test
emmeans.tsoil.19 <- emmeans(lme.tsoil.19, c("Treatment"), type = "response")
pairs(emmeans.tsoil.19, type = "response", adjust = "bonf")

# lme treatment effects on VWC with Block as a random factor on growing season averages

lme.vwc.18 <-lme(VWC_mean ~ Treatment,
              random = ~ 1|plot/Block,
              data = na.omit(subset(env_data_gs_18_19, year %in% c("18"))))
anova.vwc18 <- anova(lme.vwc.18)
plot(resid(lme.vwc.18))
qqnorm(resid(lme.vwc.18))

emmeans.vwc.18 <- emmeans(lme.vwc.18, c("Treatment"), type = "response")
pairs(emmeans.vwc.18, type = "response", adjust = "bonf")


lme.vwc.19 <-lme(VWC_mean ~ Treatment,
                 random = ~ 1|plot/Block,
                 data = na.omit(subset(env_data_gs_18_19, year %in% c("19"))))
anova.vwc19 <- anova(lme.vwc.19)
plot(resid(lme.vwc.19))
qqnorm(resid(lme.vwc.19))

#post hoc test
emmeans.vwc.19 <- emmeans(lme.vwc.19, c("Treatment"), type = "response")
pairs(emmeans.vwc.19, type = "response", adjust = "bonf")

