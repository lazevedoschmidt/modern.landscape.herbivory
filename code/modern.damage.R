#R version 3.6.1.

#load packages
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(agricolae) #You need this for the Tukey test
library(EnvStats) #rosnerTest

#color palette for ecosystems
coleco <- c("#ff6b35","#004e89","#2c6e49")

#load data file----
dam_data <- read.csv("Landsape.spreadsheet_quarry.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
dam_data$Forest <- c("La Selva","La Selva","La Selva","La Selva","La Selva",
                     "La Selva","La Selva","Harvard Forest",
                     "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                     "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                     "Harvard Forest","SERC","SERC","SERC","SERC","SERC","SERC",
                     "SERC","SERC","SERC")
#Total percent damage
avganth.all <- dam_data %>%
  #drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Forest)%>% #grouping by epoch
  summarise(mean.all= mean(perc.dam), #average
            sd.all = sd(perc.dam), #standard deviation
            mean.spec=mean(perc.spec),
            sd.spec=sd(perc.spec),
            mean.spec=mean(perc.mine),
            sd.spec=sd(perc.mine),
            mean.gall=mean(perc.gall),
            sd.gall=sd(perc.gall),
            mean.HF=mean(perc.HF),
            sd.HF=sd(perc.HF),
            mean.mf=mean(perc.mf),
            sd.mf=sd(perc.mf),
            mean.skel=mean(perc.skel),
            sd.skel=sd(perc.skel),
            mean.sf=mean(perc.sf),
            sd.sf=sd(perc.sf),
            mean.pier=mean(perc.pierce),
            sd.pier=sd(perc.pierce),
            mean.rawdt=mean(raw.dt),
            sd.rawdt=sd(raw.dt),
            mean.rawffg=mean(raw.ffg),
            sd.rawffg=sd(raw.ffg),
            mean.divdt=mean(div.dt),
            sd.divdt=sd(div.dt),
            mean.divspec=mean(div.spec),
            sd.divspec=sd(div.spec),
            mean.divmine=mean(div.mine),
            sd.divmine=sd(div.mine),
            mean.divgall=mean(div.gall),
            sd.divgall=sd(div.gall),
            mean.shannon=mean(Shannon_plant),
            sd.shannon=sd(Shannon_plant),
            mean.pj=mean(Pj_plant),
            sd.pj=sd(Pj_plant),
            mean.shannondt=mean(DT_shan),
            sd.shannont=sd(DT_shan),
            mean.pjdt=mean(DT_Pj),
            sd.pjdt=sd(DT_Pj),
            mean.divplant=mean(Leaf.diversity_rare),
            sd.divplant=sd(Leaf.diversity_rare),
            count=n())

#plotting damage percent
avganth.all$Forest <- factor(avganth.all$Forest, 
                           levels = c("La Selva", "SERC", "Harvard Forest"))
avganth.all <- arrange(avganth.all, Forest)

dam_data$Forest <- factor(dam_data$Forest, 
                             levels = c("La Selva", "SERC", "Harvard Forest"))
dam_data <- arrange(dam_data, Forest)

#Total percent damage----
total.lm <- lm(perc.dam ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av) #checking additional TUkey test 

#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.dam,
                            k = 1)
ros.totalperc 

mod.total.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.dam, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","a"),
                y = c(60,60,60), size=2))+
  ylab("Total Damage Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.total.perc

#Total damage diversity----
total.lm <- lm(div.dt ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#checking outliers:
ros.totalperc <- rosnerTest(dam_data$div.dt,
                            k = 1)
ros.totalperc

#plotting damage diversity
mod.total.div <- ggplot(dam_data, aes(fill=Forest, y=div.dt, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","a"),
                y = c(22,22,22), size=2))+
  ylab("Total Damage Diversity \n(Standardized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.total.div

#Specialization percent----
total.lm <- lm(perc.spec ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.spec,
                            k = 1)
ros.totalperc

mod.perc.spec <- ggplot(dam_data, aes(fill=Forest, y=perc.spec, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a*","b","b"),
                y = c(17,17,17), size=2))+
  ylab("Specialized Damage Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.perc.spec

#speialized diveristy----
total.lm <- lm(div.spec ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#checking outliers:
ros.totalperc <- rosnerTest(dam_data$div.spec,
                            k = 1)
ros.totalperc

#plotting damage percent
mod.spec.div <- ggplot(dam_data, aes(fill=Forest, y=div.spec, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","a"),
                y = c(7,7,7), size=2))+
  ylab("Specialized Damage Diversity \n(Standardized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.spec.div

#Minning percent----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.mine,
                            k = 1)
ros.totalperc

total.lm <- lm(perc.mine ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting damage percent
mod.perc.mine <- ggplot(dam_data, aes(fill=Forest, y=perc.mine, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","b","b"),
                y = c(-0.5,-0.5,-0.5), size=2))+
  ylab("Mining Damage Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.perc.mine

#Mining diversity----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$div.mine,
                            k = 1)
ros.totalperc

total.lm <- lm(div.mine ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting damage percent
mod.div.mine <- ggplot(dam_data, aes(fill=Forest, y=div.mine, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","ab","b"),
                y = c(-0.5,-0.5,-0.5), size=2))+
  ylab("Mine Damage Diversity \n(Standardized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.div.mine

#Gall percent----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.gall,
                            k = 1)
ros.totalperc

total.lm <- lm(perc.gall ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting damage percent
mod.gall.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.gall, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","c","b"),
                y = c(0,0,0), size=2))+
  ylab("Gall Damage Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.gall.perc

#Gall diveristy----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$div.gall,
                            k = 1)
ros.totalperc

total.lm <- lm(div.gall ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
mod.gall.div <- ggplot(dam_data, aes(fill=Forest, y=div.gall, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","b","b"),
                y = c(0.5,0.5,0.5), size=2))+
  ylab("Gall Damage Diversity \n(Standardized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.gall.div

#Hole feeding----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.HF,
                            k = 1)
ros.totalperc

total.lm <- lm(perc.HF ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
mod.HF.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.HF, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","b","a"),
                y = c(29,29,29), size=2))+
  ylab("Hole Feeding Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.HF.perc

#Margin feeding----
#checking outliers:
ros.totalperc <- rosnerTest(dam_data$perc.mf,
                            k = 1)
ros.totalperc
dam_data <- dam_data[-26,] #remove outlier but then rerun dam_data for the rest of the script

total.lm <- lm(perc.mf ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
avganth.all$Forest <- factor(avganth.all$Forest, 
                             levels = c("La Selva", "SERC", "Harvard Forest"))
avganth.all <- arrange(avganth.all, Forest)

dam_data$Forest <- factor(dam_data$Forest, 
                          levels = c("La Selva", "SERC", "Harvard Forest"))
dam_data <- arrange(dam_data, Forest)
mod.mf.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.mf, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("b","ab","a"),
                y = c(32,32,32), size=2))+
  ylab("Margin Feeding Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.mf.perc

#Skeletonization----
ros.totalperc <- rosnerTest(dam_data$perc.skel,
                            k = 1)
ros.totalperc

total.lm <- lm(perc.skel ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
mod.skel.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.skel, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("b","a","a"),
                y = c(0,0,0), size=2))+
  ylab("Skeletonized Feeding \nFrequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.skel.perc

#Surface Feeding----
ros.totalperc <- rosnerTest(dam_data$perc.sf,
                            k = 1)
ros.totalperc

total.lm <- lm(perc.sf ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting damage percent
mod.sf.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.sf, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","a"),
                y = c(10,10,10), size=2))+
  ylab("Surface Feeding Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.sf.perc

#Piercing----
ros.totalperc <- rosnerTest(dam_data$perc.pierce,
                            k = 1)
ros.totalperc
dam_data <- dam_data[-22,]
dam_data <- dam_data[-20,]
dam_data <- dam_data[-20,]
dam_data <- dam_data[-8,]
#rerun dam data orignial code

total.lm <- lm(perc.pierce ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting damage percent
avganth.all$Forest <- factor(avganth.all$Forest, 
                             levels = c("La Selva", "SERC", "Harvard Forest"))
avganth.all <- arrange(avganth.all, Forest)

dam_data$Forest <- factor(dam_data$Forest, 
                          levels = c("La Selva", "SERC", "Harvard Forest"))
dam_data <- arrange(dam_data, Forest)

mod.p.perc <- ggplot(dam_data, aes(fill=Forest, y=perc.pierce, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("b","a","ab"),
                y = c(-0.15,-0.15,-0.15), size=2))+
  ylab("Piercing and Sucking \nFeeding Frequency (%)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.p.perc

#Shannon Diversity Plant----
ros.totalperc <- rosnerTest(dam_data$Shannon_plant,
                            k = 1)
ros.totalperc

total.lm <- lm(Shannon_plant ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
avganth.all$Forest <- factor(avganth.all$Forest, 
                             levels = c("La Selva", "SERC", "Harvard Forest"))
avganth.all <- arrange(avganth.all, Forest)

dam_data$Forest <- factor(dam_data$Forest, 
                          levels = c("La Selva", "SERC", "Harvard Forest"))
dam_data <- arrange(dam_data, Forest)

mod.shannon <- ggplot(dam_data, aes(fill=Forest, y=Shannon_plant, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","b"),
                y = c(0.5,0.5,0.5), size=2))+
  ylab("Shannon Diversity (Leaf)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.shannon

#Shannon Diversity DT----
ros.totalperc <- rosnerTest(dam_data$DT_shan,
                            k = 1)
ros.totalperc

total.lm <- lm(DT_shan ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
mod.shannondt <- ggplot(dam_data, aes(fill=Forest, y=DT_shan, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","b","c"),
                y = c(0.5,0.5,0.5), size=2))+
  ylab("Shannon Diversity (DT)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.shannondt

#PJ (leaf)----
ros.totalperc <- rosnerTest(dam_data$Pj_plant,
                            k = 1)
ros.totalperc

total.lm <- lm(Pj_plant ~ Forest, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Forest')
total_test
TukeyHSD(total.av)

#plotting
mod.PJ <- ggplot(dam_data, aes(fill=Forest, y=Pj_plant, x=Forest)) +
  geom_boxplot()+
  geom_jitter(aes(shape=Dep.env), color = "black",size=4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  geom_text(data= data.frame(avganth.all), #lumped modern
            aes(label = c("a","a","a"),
                y = c(0.4,0.4,0.4), size=2))+
  ylab("Pielou's J (Evenness; Leaf)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
mod.PJ

#PJ (DT)----
# ros.totalperc <- rosnerTest(dam_data$DT_Pj,
#                             k = 1)
# ros.totalperc
# 
# total.lm <- lm(DT_Pj ~ Forest, data = dam_data)
# total.av <- aov(total.lm)
# summary(total.av)
# total_test <- HSD.test(total.av, 'Forest')
# total_test
# TukeyHSD(total.av)
# 
# #plotting
# mod.PJ.dt <- ggplot(dam_data, aes(fill=Forest, y=DT_Pj, x=Forest)) +
#   geom_boxplot()+
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
#   scale_fill_manual(values=coleco)+ #color palette created above
#   #ggtitle("Average Percent Damage by Epoch")+ #title
#   geom_text(data= data.frame(avganth.all), #lumped modern
#             aes(label = c("a","a","a"),
#                 y = c(0.7,0.7,0.7), size=2))+
#   ylab("Pielou's J (Evenness; DT)")+ #axis labels
#   xlab("")+
#   theme_light()+ #light background makes bar plot easier to read
#   theme(legend.position = "none") #removes legend
# mod.PJ.dt

#making total super plot to show no difference across all eosystems----
#Shannon diversity and evenenness: plant
super.shan <- plot_grid(mod.shannon,
                        mod.PJ,
                        nrow = 1,
                        ncol = 2,
                        labels = c("A", "B"))
super.shan

super.total_spec <- plot_grid(mod.total.div, mod.total.perc,
                              mod.spec.div, mod.perc.spec,
                         nrow = 1,
                         ncol = 4, 
                         labels = c("C", "D","E", "F"))
super.total_spec

#specialist FFG figure
mine_gall <- plot_grid(mod.perc.mine,
                      mod.div.mine,
                      mod.gall.perc,
                      mod.gall.div,
                      nrow = 1,
                      ncol = 4, 
                      labels = c("G", "H", "I", "J"))
mine_gall

spec.sig <- plot_grid(mod.sf.perc,
                      mod.skel.perc,
                      mod.p.perc,
                      nrow = 1, 
                      ncol = 3, 
                      labels = c("K", "L","M"))
spec.sig

#gen. sig figure
super.gen.sig <- plot_grid(mod.HF.perc,
                           mod.mf.perc,
                           nrow=1,
                          ncol = 2,
                          labels = c("N", "O"))
super.gen.sig


#making full figure for manuscript
fullmodplot <- plot_grid(super.shan,
                         super.total_spec,
                         mine_gall,
                         spec.sig,
                         super.gen.sig,
                         nrow = 5,
                         ncol = 1)
fullmodplot
#Saving figures----
ggsave("figures/fullmodplot.pdf", fullmodplot, width = 10, height = 15, units = "in")
