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
boot_data <- read.csv("bootstrapmean.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
boot_data$Forest <- c("La Selva","La Selva","La Selva","La Selva","La Selva",
            "La Selva","La Selva","Harvard Forest",
            "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
            "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
            "Harvard Forest","SERC","SERC","SERC","SERC","SERC","SERC",
            "SERC","SERC","SERC")
#ordering correctly for plotting 
boot_data$Forest <- factor(boot_data$Forest, 
                           levels = c("La Selva", "SERC", "Harvard Forest"))
boot_data <- arrange(boot_data, Forest)
boot_data <- boot_data %>%
  mutate(nest2 = nestedness/10) %>%
  relocate(nest2, .before="Forest")

hfls <- filter(boot_data, Forest == "La Selva" | Forest == "Harvard Forest")
hfmd <- filter(boot_data, Forest == "Harvard Forest" | Forest == "SERC")

#making dataframe long for raincloud plot
long_boot <-  gather(data=boot_data, key=ID, value=value, connectance:nest2, 
                     factor_key = TRUE)
long_boot$X <- NULL


#HF/LS comparisons
lm <- lm(connectance ~ Forest, data = hfls)
av <- aov(lm)
summary(av)
total_test <- HSD.test(av, 'Forest')
total_test
TukeyHSD(av) #checking additional TUkey test 

lm <- lm(robustness.LL ~ Forest, data = hfls)
av <- aov(lm)
summary(av)
total_test <- HSD.test(av, 'Forest')
total_test
TukeyHSD(av)

connect <- ggplot(boot_data, aes(fill=Forest, y=connectance, x=Forest)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  scale_fill_manual(values=coleco)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  # geom_text(data= data.frame(avganth.all), #lumped modern
  #           aes(label = c("a","a","a"),
  #               y = c(60,60,60), size=2))+
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
connect 

#updated: 03/07/2022
#raincloud plots for HF-LS
col2 <- c("#ff6b35","#2c6e49")
HF.LS.rain <- long_boot %>%
  filter(Forest == "La Selva" | Forest == "Harvard Forest") %>%
  filter(ID == "connectance"| ID == "H2" | ID == "partner.diversity.HL" | ID == "nest2") %>%
  ggplot(aes(x = value, y = ID, fill= Forest,
             color = Forest))+
  ggdist::stat_halfeye(adjust = 1, # this changes the intervals used for calculating the density plot (bigger=smoother)
                       justification = -0.25,
                       .width = 0,
                       width = 0.8, 
                       point_colour = NA,
                       trim=FALSE, # whether to trim to data range-- looks nicer without trimming, but can be misleading
                       normalize="xy", # this normalizes each line, rather than having them all on the same scale
                       scale=0.6) + # change to adjust height (scale=1 means that they fill the whole line height)
  geom_point(
    size = 1.3,
    alpha = 0.3,
    position = position_jitter(
      seed = 1, width=0.001, height = 0.1
    )
  ) +
  geom_boxplot(width = 0.3,
               outlier.colour = NA, 
               alpha = 0.25,
               position = position_dodge(width=0)) + # this makes the boxplots overlap-- remove if you'd rather they be side-by-side :)
  # ggdist::stat_dots(side= "left",
  #                   justification = 1.15,
  #                   binwidth = 0.5)+
  
  #coord_flip()+
  scale_fill_manual(values=col2)+
  scale_color_manual(values=col2)+
  theme_light()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16)) +
  #legend.position = c(0.8,0.9), #add these back in if we want to see the legend 
  #legend.title = element_blank(),
  #legend.text = element_text(size=16),
  #coord_cartesian(xlim=c(0,2.5),clip="off") +
  #scale_y_discrete(labels= freqlabels)+
  theme(legend.position = "none") +
  #theme(legend.position = c(0.87, 0.25))+
  xlab("")
HF.LS.rain

#nestedness needs to be alone due to scale issue
# HF.LS.rain.nest <- long_boot %>%
#   filter(Forest == "La Selva" | Forest == "Harvard Forest") %>%
#   filter(ID == "nestedness") %>%
#   ggplot(aes(x = value, y = ID, fill= Forest,
#              color = Forest))+
#   ggdist::stat_halfeye(adjust = 1, # this changes the intervals used for calculating the density plot (bigger=smoother)
#                        justification = -0.25,
#                        .width = 0,
#                        width = 0.8, 
#                        point_colour = NA,
#                        trim=FALSE, # whether to trim to data range-- looks nicer without trimming, but can be misleading
#                        normalize="xy", # this normalizes each line, rather than having them all on the same scale
#                        scale=0.6) + # change to adjust height (scale=1 means that they fill the whole line height)
#   geom_point(
#     size = 1.3,
#     alpha = 0.3,
#     position = position_jitter(
#       seed = 1, width=0.001, height = 0.1
#     )
#   ) +
#   geom_boxplot(width = 0.3,
#                outlier.colour = NA, 
#                alpha = 0.25,
#                position = position_dodge(width=0)) + # this makes the boxplots overlap-- remove if you'd rather they be side-by-side :)
#   # ggdist::stat_dots(side= "left",
#   #                   justification = 1.15,
#   #                   binwidth = 0.5)+
#   
#   #coord_flip()+
#   scale_fill_manual(values=col2)+
#   scale_color_manual(values=col2)+
#   theme_light()+
#   theme(axis.text = element_text(size=14),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size=16)) +
#   #legend.position = c(0.8,0.9), #add these back in if we want to see the legend 
#   #legend.title = element_blank(),
#   #legend.text = element_text(size=16),
#   #coord_cartesian(xlim=c(0,2.5),clip="off") +
#   #scale_y_discrete(labels= freqlabels)+
#   #theme(legend.position = "none") +
#   theme(legend.position = c(0.87, 0.25))+
#   xlab("")
# HF.LS.rain.nest

#SERC-HF
cole2v2 <- c("#004e89","#2c6e49")
HF.MD.rain <- long_boot %>%
  filter(Forest == "SERC" | Forest == "Harvard Forest") %>%
  filter(ID == "connectance"| ID == "H2" | ID == "partner.diversity.HL"| ID == "robustness.LL") %>%
  ggplot(aes(x = value, y = ID, fill= Forest,
             color = Forest))+
  ggdist::stat_halfeye(adjust = 1, # this changes the intervals used for calculating the density plot (bigger=smoother)
                       justification = -0.25,
                       .width = 0,
                       width = 0.8, 
                       point_colour = NA,
                       trim=FALSE, # whether to trim to data range-- looks nicer without trimming, but can be misleading
                       normalize="xy", # this normalizes each line, rather than having them all on the same scale
                       scale=0.6) + # change to adjust height (scale=1 means that they fill the whole line height)
  geom_point(
    size = 1.3,
    alpha = 0.3,
    position = position_jitter(
      seed = 1, width=0.001, height = 0.1
    )
  ) +
  geom_boxplot(width = 0.3,
               outlier.colour = NA, 
               alpha = 0.25,
               position = position_dodge(width=0)) + # this makes the boxplots overlap-- remove if you'd rather they be side-by-side :)
  # ggdist::stat_dots(side= "left",
  #                   justification = 1.15,
  #                   binwidth = 0.5)+
  
  #coord_flip()+
  scale_fill_manual(values=cole2v2)+
  scale_color_manual(values=cole2v2)+
  theme_light()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16)) +
  #legend.position = c(0.8,0.9), #add these back in if we want to see the legend 
  #legend.title = element_blank(),
  #legend.text = element_text(size=16),
  #coord_cartesian(xlim=c(0,2.5),clip="off") +
  #scale_y_discrete(labels= freqlabels)+
  #theme(legend.position = "none") +
  theme(legend.position = c(0.87, 0.25))+
  xlab("")
HF.MD.rain

#making plots
# full.hf_ls.plot <- plot_grid(HF.LS.rain, 
#                         HF.LS.rain.nest,
#                         nrow = 1, 
#                         ncol = 2, 
#                         labels = c("A", "B"))
# full.hf_ls.plot

all.plots <- plot_grid(HF.LS.rain,
                       HF.MD.rain,
                       nrow = 2,
                       ncol = 1,
                       labels = c("A", "B"))
all.plots

#saving figures: 
ggsave("figures/full_rain.pdf", all.plots, width = 12.5, height = 15,units = "in")



