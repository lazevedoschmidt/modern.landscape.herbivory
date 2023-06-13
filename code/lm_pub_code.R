#R version 3.6.1.
#Models for publication

#load packages: 
library(lme4) #lmer 
library(ggplot2)
library(cowplot)
library(car) #needed for logit transformation
library(languageR) #needed for p values of models using pvals.fnc(model) function
library(tidyverse)
library(effects) #partial residual plots
library(EnvStats)

#needed to run sjPlot
library(insight)
library(performance) 
library(parameters)
library(effectsize)
library(sjPlot) #needed for making lmer summary tables
library(webshot) #needed to save sjPlot tables

dam_data <- read.csv("Landsape.spreadsheet_quarry.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
dam_data <- dam_data %>%
  mutate(Forest = case_when(Environment == "temperate" ~ "Harvard Forest",
                            Environment == "coastal temp" ~ "SERC", 
                            Environment == "tropical wet" ~ "La Selva"))
dam_data$note <- NULL
dam_data$Lat <- scale(dam_data$Lat, center = T, scale = T) 
dam_data$Elevation..m. <- scale(dam_data$Elevation..m., center = T, scale = T)
dam_data$Shannon_plant <- scale(dam_data$Shannon_plant, center = T, scale = T)
dam_data$Pj_plant <- scale(dam_data$Pj_plant, center = T, scale = T)

dam_data$perc.dam_logit <- logit(dam_data$perc.dam, percents = TRUE)
dam_data$perc.gall_logit <- logit(dam_data$perc.gall, percents = TRUE) 
dam_data$perc.HF_logit <- logit(dam_data$perc.HF, percents = TRUE)
dam_data$perc.mf_logit <- logit(dam_data$perc.mf, percents = TRUE)
dam_data$perc.mine_logit <- logit(dam_data$perc.mine, percents = TRUE) 
#what does this warning mean: 
#In logit(dam_data$perc.mine, percents = TRUE) :
#proportions remapped to (0.025, 0.975)
dam_data$perc.pierce_logit <- logit(dam_data$perc.pierce, percents = TRUE) 
dam_data$perc.sf_logit <- logit(dam_data$perc.sf, percents = TRUE)
dam_data$perc.skel_logit <- logit(dam_data$perc.skel, percents = TRUE)
dam_data$perc.spec_logit <- logit(dam_data$perc.spec, percents = TRUE)

mod1 <- lmer(perc.dam_logit ~ Lat + Shannon_plant + Pj_plant +
               (1|Forest/Facies) , data=dam_data) 
summary(mod1)
AIC(mod1)

mod1.1 <- lmer(perc.spec_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.1)
AIC(mod1.1)

mod1.2 <- lmer(perc.gall_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.2)
AIC(mod1.2)

mod1.3 <- lmer(perc.mine_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.3)
AIC(mod1.3)

mod1.4 <- lmer(perc.HF_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.4)
AIC(mod1.4)

mod1.5 <- lmer(perc.mf_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.5)
AIC(mod1.5)

mod1.6 <- lmer(perc.skel_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.6)
AIC(mod1.6)

mod1.7 <- lmer(perc.pierce_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.7)
AIC(mod1.7)

mod1.8 <- lmer(perc.sf_logit ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.8)
AIC(mod1.8)

mod1.9 <- lmer(div.dt ~ Lat + Shannon_plant + Pj_plant +
                 (1|Forest/Facies) , data=dam_data) 
summary(mod1.9)
AIC(mod1.9)

mod1.10 <- lmer(div.spec ~ Lat + Shannon_plant + Pj_plant +
                  (1|Forest/Facies) , data=dam_data) 
summary(mod1.10)
AIC(mod1.10)

mod1.11 <- lmer(div.gall ~ Lat + Shannon_plant + Pj_plant +
                  (1|Forest/Facies) , data=dam_data) 
summary(mod1.11)
AIC(mod1.11)

mod1.12 <- lmer(div.mine ~ Lat + Shannon_plant + Pj_plant +
                  (1|Forest/Facies) , data=dam_data) 
summary(mod1.12)
AIC(mod1.12)

#Summary tables: model 1----
tab_model(mod1,mod1.1,mod1.2,mod1.3,
          show.re.var= TRUE,
          show.aic = TRUE,
          pred.labels =c("(Intercept)", "Latitude","Plant Diversity (Shannon)", "Plant Evenness"),
          #dv.labels= "Harvard Forest",
          file = "tables/Latlmer.freq.sumtab.pt1.html")

tab_model(mod1.4,mod1.5,mod1.6,mod1.7,mod1.8,
          show.re.var= TRUE,
          show.aic = TRUE,
          pred.labels =c("(Intercept)", "Latitude","Plant Diversity (Shannon)", "Plant Evenness"),
          #dv.labels= "Harvard Forest",
          file = "tables/Latlmer.freq.sumtab.pt2.html")


tab_model(mod1.9,mod1.10,mod1.11,mod1.12,
          show.re.var= TRUE,
          show.aic = TRUE,
          pred.labels =c("(Intercept)", "Latitude","Plant Diversity (Shannon)", "Plant Evenness"),
          #dv.labels= "Harvard Forest",
          file = "tables/Latlmer.div.sumtab.html")