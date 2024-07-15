# "Jädraås Respiration Paper Figure 1a-c"
#May 3 2021
#Updated 5 May 2021
rm(list=ls()) #clear the workspace


#libraries
library(mgcv); library(nlme); library(ggplot2); library(multcomp); library(car); library(ape);library(readxl);
library(itsadug);library(MASS);library(grid);library(readxl);library(ggplot2)
library(ggsci);library(tidyr);library(tidyverse);library(scales);library(gridExtra)
library(grid);library(lattice);library(MASS);library(lme4);library(nlme);library(emmeans);library(gamm4);library(dplyr);library(magrittr);library(reshape2)
library(mgcv);library(lmerTest);library(lubridate);library(car);library(lmerTest);library(pls); library(forcats);library(patchwork)


# loading all soil temp and vwc data

setwd("/Users/lske0002/Projects/mycorrhizal removal/Soil VWC & Temp/data")
env.data <- read.csv("loggers_2016Dec07-2019Nov20_long.csv")
env.data <- na.omit(env.data)

#characterizing time categories
env.data$Date <- as.Date(env.data$Date)
env.data$Treatment <- factor(env.data$Treatment, levels = c('SPE', 'SP', 'SE','S','DC'))

env.data$d <- as.numeric(decimal_date(env.data$Date)) # making a vector of years with fractions

# setting negative values of moisture to 0
env.data$VWC <- ifelse(env.data$VWC < 0, 0, env.data$VWC) 

#making Block a factor, not an integer & setting up a column with plot ID
env.data$Block <- as.factor(env.data$Block)
env.data$plot <- as.factor(paste(env.data$Treatment,env.data$Block))  

#summary for vwc data for each day over all plots
vwc_jädraås1 <- env.data %>%
  group_by(Date, Block) %>% 
  drop_na() %>% 
  summarize_at(c("VWC"), funs(mean,n(),se=sd(.)/sqrt(n())))
vwc_jädraås1

vwc_jädraås1$Treatment <- rep("All",1079) #dummy variables
vwc_jädraås1$Year <- rep("NA",1079) #dummy variables


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

tsoil_jädraås1$Treatment <- rep("All",1079) #dummy variables
tsoil_jädraås1$Year <- rep("NA",1079) #dummy variables

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

######### Figure 1 a Respiration Point Measurements
data.gs <- read.csv("/Users/lske0002/Projects/mycorrhizal removal/Response Curves/data/Tsoil&VWC_fourdayavg_block&wholexpgapfill_noNov_13JULY2020.csv") 

#growing season (no Nov data)
data.gs$Block <- as.factor(data.gs$Block) 
data.gs$erm_shrubs <- as.factor(data.gs$erm_shrubs)
data.gs$ecm_roots <- as.factor(data.gs$ecm_roots)
data.gs$Date <- as.Date(data.gs$Date)
data.gs$Month <- as.numeric(format(data.gs$Date, format = "%m"))
data.gs$d <- decimal_date(data.gs$Date) #making a vector of years with fractions
data.gs$Year1<- ifelse(data.gs$Year=="2017",1,ifelse(data.gs$Year=="2018",2,3))
data.gs$plot <- paste(data.gs$Treatment,data.gs$Block, sep = "")

#subset by treatment
data.gs_DC <- data.gs[data.gs$Treatment == "DC",]
data.gs_S <- data.gs[data.gs$Treatment == "S",] # n = 200
data.gs_SE <- data.gs[data.gs$Treatment == "SE",] # n = 200
data.gs_SP <- data.gs[data.gs$Treatment == "SP",] # n = 200
data.gs_SPE <- data.gs[data.gs$Treatment == "SPE",] # n = 200

#combining treatments based on downstream analyses
data.gs_noDC <- rbind(data.gs_S,data.gs_SE,data.gs_SP,data.gs_SPE)
#combining treatments for graphing
data.gs_noDC$Respiration_mg <- (data.gs_noDC$Respiration)*(12.01/44.01)*1000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# averaging soil respiration data across blocks by each date measurements were taken
data_grouped <- data.gs_noDC %>%
  group_by(Year,Date,Treatment) %>%
  drop_na() %>% 
  summarize_at(c("Respiration_mg"), funs(mean, sd, n(), se=sd(.)/sqrt(n())))
resp.gs <- ungroup(data_grouped)

# Define Start and end times for the subset as R objects that are the time class
startTime <- as.Date("2017-01-1")
endTime <- as.Date("2019-12-01")

start.end <- c(startTime,endTime)
start.end

#colors
cbp2 <- c("#D55E00","#E6CF02","#009E73", "#56B4E9", "#E69F00")

#plot growing season Jädraås mycorrhizal roots and saprotrophic respiration rates 2017-2019
df <- data.frame(grp = resp.gs$Date, fit = resp.gs$mean, se = resp.gs$se)
figure1a_resp <- ggplot(df, aes(grp, fit, ymin = fit-se, ymax = fit+se, color=resp.gs$Treatment)) +
  geom_line(aes(color = resp.gs$Treatment))+
  geom_pointrange(size = 1, shape = 1) +
  xlab(element_blank()) +
  ylab(expression(Soil*" "*respiration*" "*mg*" "*C*" "*m^-2*" "*h^-1)) +
  scale_x_date(limits=start.end, date_breaks = "2 months" , date_labels = "%b")+
  scale_y_continuous(expand =c(0,0),limits = c(0,415))+
  theme(axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        #legend.title = element_text(size = 12),
        #legend.text = element_text(size = 12),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        panel.border= element_rect(colour = "black", size=1, fill=NA),
        panel.background = element_rect(colour = "black", size=.5, fill=NA)) +
  scale_shape_manual(values=c(1,1,1,1))+
  scale_color_manual(values = c(cbp2))+
  guides(color=guide_legend(title="Treatment"))
figure1a_resp


##############-------------------------------##########

############## soil moist point measurements ##########

# averaging soil moisture by block
data.vwc.point.avg <- data.gs_noDC %>%
  group_by(Treatment, Year, Date) %>%
  drop_na() %>% 
  summarize_at(c("VWC_mean"), funs(mean, n(), se=sd(.)/sqrt(n())))
data.vwc.point.avg <- ungroup(data.vwc.point.avg)
data.vwc.point.avg <- data.vwc.point.avg[-1,]

#adding the all plot average of soil moisture to fill out the plot
alla_moist <- rbind(data.vwc.point.avg,vwc_jädraås1)

#colors
cbp <- c("black","#E6CF02", "#009E73","#56B4E9",
         "#E69F00")
#plot 
df <- data.frame(grp = alla_moist$Date, fit = alla_moist$mean, se = alla_moist$se, treatment = alla_moist$Treatment)
figure1b_moist <- ggplot(df, aes(shape = "1", grp, fit, ymin = fit-se, ymax = fit+se, color=treatment)) +
  theme_classic() +
  ggtitle("     ") +
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
  geom_line(data = df[df$treatment == "All",])+
  geom_pointrange(data = df[df$treatment != "All",],size = 1) +
  scale_color_manual(values = c(cbp))+
  scale_shape_manual(values=c(1,1,1,1))+
scale_fill_discrete(guide = guide_legend(reverse=TRUE))
figure1b_moist

############# soil temp point measurements ###############

# averaging soil temp across blocks within treatments, dates and years
data.tsoil.point.avg <- data.gs_noDC %>%
  group_by(Treatment, Year, Date) %>%
  drop_na() %>% 
  summarize_at(c("Tsoil_mean"), funs(mean, n(), se=sd(.)/sqrt(n())))
data.tsoil.point.avg <- ungroup(data.tsoil.point.avg)
data.tsoil.point.avg <- data.tsoil.point.avg[-1,]

#adding the all plot average of soil temp to fill out the plot
alla_tsoil <- rbind(data.tsoil.point.avg,tsoil_jädraås1)

#plot 
df <- data.frame(grp = alla_tsoil$Date, fit = alla_tsoil$mean, se = alla_tsoil$se, treatment = alla_tsoil$Treatment)
figure1c_temp <- ggplot(df, aes(shape = "1", grp, fit, ymin = fit-se, ymax = fit+se, color=treatment)) +
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
figure1c_temp

#combining plots into display
require(gridExtra)
grid.arrange(figure1a_resp, figure1b_moist, figure1c_temp, nrow=3, padding = c(0,0))



#add month as categorical variable
env.data$month <- as.factor(as.numeric((format(env.data$Date, format = "%m"))))
#subset the dataset by growing season months (may-october)
climate.jädraås.growingseason <- subset(env.data, month %in% c("5","6","7","8","9","10"))

#add year as a categorical variable
climate.jädraås.growingseason$year <- as.factor(as.numeric((format(climate.jädraås.growingseason$Date, format = "%y"))))

#subset the dataset by growing season months (may-october)
climate.jädraås.growingseason <- subset(climate.jädraås.growingseason, year %in% c("17","18","19"))

#subset the dataset by growing season months (may-october)
climate.jädraås.growingseason <- subset(climate.jädraås.growingseason, Block %in% c(2,5,7))

#summary of jädraås climate data by year, month, Treatment and Block from daily averages
climate_jädraås_growingseason <- climate.jädraås.growingseason %>%
  group_by(Block, year) %>%
  summarize_at(c("VWC","Tsoil"), funs(mean, min, max, sd, n()))

######


