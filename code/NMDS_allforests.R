#NMDS plots by quarry for each forest

require(ggplot2)
require(vegan)
require(cowplot)
require(tidyverse)
require(numDeriv)
require(ggrepel)
require(MASS)

#load DT data file----
dam_data <- read.csv("Landsape.spreadsheet_quarry.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
dam_data <- dam_data %>%
  mutate(Forest = case_when(Environment == "temperate" ~ "Harvard Forest",
                            Environment == "coastal temp" ~ "SERC", 
                            Environment == "tropical wet" ~ "La Selva"))
dam_data$note <- NULL

dam_data2 <- dam_data %>%
  filter(Forest == "Harvard Forest" | Forest == "SERC")

dam_data3 <- dam_data %>%
  filter(Forest == "La Selva")

#creating proportion of each FFG damage
#HF and SERc only----
propdam <- dam_data2 %>%
  group_by(Facies)%>%
  summarise(summine = sum(perc.mine), 
            sumgall = sum(perc.gall),
            sumHF = sum(perc.HF),
            summf = sum(perc.mf),
            sumskel = sum(perc.skel),
            sumsf = sum(perc.sf),
            sumpierce = sum(perc.pierce))

newdam <- full_join(dam_data2, propdam, by = "Facies")

newdam2 <- newdam %>%
  mutate(mine = perc.mine/summine,
         gall = perc.gall/sumgall,
         HF = perc.HF/sumHF,
         MF = perc.mf/summf,
         skel = perc.skel/sumskel,
         SF = perc.sf/sumsf,
         P.S = perc.pierce/sumpierce)

#La Selva only: without p/s damage----
propdam2 <- dam_data3 %>%
  group_by(Facies)%>%
  summarise(summine = sum(perc.mine), 
            sumgall = sum(perc.gall),
            sumHF = sum(perc.HF),
            summf = sum(perc.mf),
            sumskel = sum(perc.skel),
            sumsf = sum(perc.sf))

newdam3 <- full_join(dam_data3, propdam2, by = "Facies")

newdam3 <- newdam3 %>%
  mutate(mine = perc.mine/summine,
         gall = perc.gall/sumgall,
         HF = perc.HF/sumHF,
         MF = perc.mf/summf,
         skel = perc.skel/sumskel,
         SF = perc.sf/sumsf)

#Cleaning up data for NMDS (vegan) analyses
newdam_hfmd <- newdam2 %>%
  dplyr::select(Quarry, Forest, Dep.env,mine,gall,HF,MF,
                skel,SF,P.S) 
#total,specialized,divtot,divspec,divgall,
#divmine,propplant

newdam_ls <- newdam3 %>%
  dplyr::select(Quarry, Forest, Dep.env, mine,gall,HF,MF,
                skel,SF)
#total,specialized,divtot,divspec,divgall,
#divmine,propplant

allfor <- full_join(newdam_hfmd,newdam_ls)
allfor$Forest <- NULL
allfor$Dep.env <- NULL
rownames(allfor) <- allfor$Quarry
allfor$Quarry <- NULL
allfor[is.na(allfor)] = 0 #making LS P/s damage NAs to 0


#DT Analysis: ALL forests----
sol.full <- metaMDS(allfor, distance = "bray", k=2, try = 20, trymax = 200)
sol.full
plot(sol.full, type="t") #Done with 2003 and 2014 data in dataframe
stressplot(sol.full)

#NMDS in ggplot
data.scores <- as.data.frame(scores(sol.full))  
#Using the scores function from vegan to extract the site scores 
#and convert to a data.frame
data.scores$Quarry <- rownames(data.scores) 
# create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(sol.full, "species"))  
#Using the scores function from vegan to extract the species scores 
#and convert to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
head(species.scores)

data.scores$dep.env <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                         "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland",
                         "Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                         "Tributary", "Fluvial", "Fluvial", "Fluvial","Fluvial", "Fluvial", "Fluvial", "Tributary", "Swamp", 
                         "Swamp", "Swamp")
data.scores$forest <- c("Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                        "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                        "SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC",
                        "La Selva","La Selva","La Selva","La Selva","La Selva","La Selva","La Selva")

jitter <- position_jitter(width = 0.1, height = 0.1)

NMDS.dt.full <- ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2)) +
  geom_text_repel(data= species.scores,aes(label = species),
                  segment.size = 0.2,
                  point.padding = NA,
                  box.padding = 0.25,
                  nudge_y = 0.1) +
  geom_point(data=species.scores, aes(x=NMDS1,y=NMDS2), shape =4)+
  geom_point(aes(color=forest), shape = c(17,17,17,15,15,15,16,16,16,3,
                                          17,17,17,15,15,15,16,16,16,
                                          16,16,16,15,17,17,17),
             size=4, position = jitter) +
  scale_color_manual(values = c("#2c6e49","#ff6b35","#004e89"))+
  #guides(colour=guide_legend(reverse = TRUE))+
  #coord_equal() +
  theme_bw() +
  #ggtitle("Insect Communities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12))+
  theme(legend.position = "none")
NMDS.dt.full

fulllegend <- get_legend(NMDS.dt.full)

facies <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland",
            "Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial","Fluvial", "Fluvial", "Fluvial", "Tributary", "Swamp", 
            "Swamp", "Swamp")
facies.hf <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland")
facies.md <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial")
facies.ls <- c("Fluvial", "Fluvial", "Fluvial", "Tributary", "Swamp", 
            "Swamp", "Swamp")
full.dt.anosim <- anosim(allfor, facies, distance = "bray") #all forests
summary(full.dt.anosim)
hf.dt.anosim <- anosim(allfor[c(1:10),], facies.hf, distance = "bray") #HF only
summary(hf.dt.anosim)
md.dt.anosim <- anosim(allfor[c(11:19),], facies.md, distance = "bray") #MD only
summary(md.dt.anosim)
ls.dt.anosim <- anosim(allfor[c(20:26),], facies.ls, distance = "bray")
summary(ls.dt.anosim)

forest <- c("Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
            "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
            "SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC",
            "La Selva","La Selva","La Selva","La Selva","La Selva","La Selva","La Selva")
full.dt.anosim2 <- anosim(allfor, forest, distance = "bray")
summary(full.dt.anosim2)

facies <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland",
            "Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
            "Tributary", "Fluvial", "Fluvial", "Fluvial",
            "Fluvial", "Fluvial", "Fluvial", "Tributary", "Swamp", 
            "Swamp", "Swamp")
full.dt.anosim3 <- anosim(allfor, facies, distance = "bray")
summary(full.dt.anosim3)

#load PLANT data file----
plant_data.hf <- read.csv("Plant_abundance_HF.csv", header=T, sep=",", 
                          na.strings=  "NA", dec=".", strip.white=TRUE)
plant_data.md <- read.csv("Plant_abundance_MD.csv", header=T, sep=",", 
                          na.strings=  "NA", dec=".", strip.white=TRUE)
plant_data.ls <- read.csv("Plant_abundance_LS.csv", header=T, sep=",", 
                          na.strings=  "NA", dec=".", strip.white=TRUE)

propplant.hf <- plant_data.hf %>%
  group_by(Facies)%>%
  summarise(sumBC = sum(BC), 
            sumBE = sum(BE),
            sumBI = sum(BI),
            sumHF1 = sum(HF1),
            sumHF2 = sum(HF2),
            sumHF3 = sum(HF3),
            sumHF4 = sum(HF4),
            sumHF5 = sum(HF5),
            sumHF6 = sum(HF6),
            sumHF7 = sum(HF7),
            sumHF8 = sum(HF8),
            sumHF9 = sum(HF9),
            sumHF14 = sum(HF14),
            sumHF15 = sum(HF15),
            sumHF16 = sum(HF16),
            sumHF17 = sum(HF17),
            sumLTA = sum(LTA),
            sumRM = sum(RM),
            sumRO = sum(RO),
            sumSU = sum(SU),
            sumWO = sum(WO))

newplantHF <- full_join(plant_data.hf, propplant.hf, by = "Facies")

newplantHF2 <- newplantHF %>%
  mutate(BC = BC/sumBC,
         BE = BE/sumBE,
         BI = BI/sumBI,
         HF1 = HF1/sumHF1,
         HF2 = HF2/sumHF2,
         HF3 = HF3/sumHF3,
         HF4 = HF4/sumHF4,
         HF5 = HF5/sumHF5,
         HF6 = HF6/sumHF6,
         HF7 = HF7/sumHF7,
         HF8 = HF8/sumHF8,
         HF9 = HF9/sumHF9,
         HF14 = HF14/sumHF14,
         HF15 = HF15/sumHF15,
         HF16 = HF16/sumHF16,
         HF17 = HF17/sumHF17,
         LTA = LTA/sumLTA,
         RM = RM/sumRM,
         RO = RO/sumRO,
         SU = SU/sumSU,
         WO = WO/sumWO)
newplantHF2[is.na(newplantHF2)] = 0

propplant.md <- plant_data.md %>%
  group_by(Facies) %>%
  summarise(sumAE = sum(AE),
            sumAH = sum(AH),
            sumBE = sum(BE),
            sumBO = sum(BO),
            sumBW = sum(BW),
            sumCO = sum(CO),
            sumGA = sum(GA),
            sumHF3 = sum(HF3),
            sumIO = sum(IO),
            sumLB = sum(LB),
            sumLTA = sum(LTA),
            sumMD1 = sum(MD1),
            sumMD2 = sum(MD2),
            sumMD3 = sum(MD3),
            sumMD4 = sum(MD4),
            sumMD5 = sum(MD5),
            sumMD6 = sum(MD6),
            sumMD7 = sum(MD7),
            sumMD8 = sum(MD8),
            sumMD9 = sum(MD9),
            sumMD10 = sum(MD10),
            sumMD11 = sum(MD11),
            sumMD12 = sum(MD12),
            sumMD13 = sum(MD13),
            sumMD14 = sum(MD14),
            sumMD15 = sum(MD15),
            sumMD16 = sum(MD16),
            sumMD17 = sum(MD17),
            sumMD18 = sum(MD18),
            sumMD19 = sum(MD19),
            sumMD20 = sum(MD20),
            sumMD21 = sum(MD21),
            sumMD22 = sum(MD22),
            sumMD23 = sum(MD23),
            sumMH = sum(MH),
            sumPIO = sum(PIO),
            sumPO = sum(PO),
            sumQS = sum(QS),
            sumRA = sum(RA),
            sumRM = sum(RM),
            sumRO = sum(RO),
            sumSCO = sum(SCO),
            sumSG = sum(SG),
            sumSRO = sum(SRO),
            sumSW = sum(SW),
            sumTP = sum(TP),
            sumWIO = sum(WIO),
            sumWO= sum(WO))

newplantMD <- full_join(plant_data.md, propplant.md, by = "Facies")

newplantMD2 <- newplantMD %>%
  mutate(AE = AE/sumAE,
         AH = AH/sumAH,
         BE = BE/sumBE,
         BO = BO/sumBO,
         BW = BW/sumBW,
         CO = CO/sumCO,
         GA = GA/sumGA,
         HF3 = HF3/sumHF3,
         IO = IO/sumIO,
         LB = LB/sumLB,
         LTA = LTA/sumLTA,
         MD1 = MD1/sumMD1,
         MD2 = MD2/sumMD2,
         MD3 = MD3/sumMD3,
         MD4 = MD4/sumMD4,
         MD5 = MD5/sumMD5,
         MD6 = MD6/sumMD6,
         MD7 = MD7/sumMD7,
         MD8 = MD8/sumMD8,
         MD9 = MD9/sumMD9,
         MD10 = MD10/sumMD10,
         MD11 = MD11/sumMD11,
         MD12 = MD12/sumMD12,
         MD13 = MD13/sumMD13,
         MD14 = MD14/sumMD14,
         MD15 = MD15/sumMD15,
         MD16 = MD16/sumMD16,
         MD17 = MD17/sumMD17,
         MD18 = MD18/sumMD18,
         MD19 = MD19/sumMD19,
         MD20 = MD20/sumMD20,
         MD21 = MD21/sumMD21,
         MD22 = MD22/sumMD22,
         MD23 = MD23/sumMD23,
         MH = MH/sumMH,
         PIO = PIO/sumPIO,
         PO = PO/sumPO,
         QS = QS/sumQS,
         RA = RA/sumRA,
         RM = RM/sumRM,
         RO = RO/sumRO,
         SCO = SCO/sumSCO,
         SG = SG/sumSG,
         SRO = SRO/sumSRO,
         SW = SW/sumSW,
         TP = TP/sumTP,
         WIO = WIO/sumWIO,
         WO= WO/sumWO)
newplantMD2[is.na(newplantMD2)] = 0

propplant.ls <- plant_data.ls %>%
  group_by(Facies)%>%
  summarise(sumLS1 = sum(LS1),
            sumLS10 = sum(LS10),
            sumLS11 = sum(LS11),
            sumLS12 = sum(LS12),
            sumLS13 = sum(LS13),
            sumLS14 = sum(LS14),
            sumLS15 = sum(LS15),
            sumLS16 = sum(LS16),
            sumLS17 = sum(LS17),
            sumLS18 = sum(LS18),
            sumLS19 = sum(LS19),
            sumLS20 = sum(LS20),
            sumLS21 = sum(LS21),
            sumLS22 = sum(LS22),
            sumLS23 = sum(LS23),
            sumLS24 = sum(LS24),
            sumLS26 = sum(LS26),
            sumLS27 = sum(LS27),
            sumLS28 = sum(LS28),
            sumLS29 = sum(LS29),
            sumLS30 = sum(LS30),
            sumLS31 = sum(LS31),
            sumLS32 = sum(LS32),
            sumLS33 = sum(LS33),
            sumLS34 = sum(LS34),
            sumLS35 = sum(LS35),
            sumLS36 = sum(LS36),
            sumLS37 = sum(LS37),
            sumLS38 = sum(LS38),
            sumLS39 = sum(LS39),
            sumLS40 = sum(LS40),
            sumLS41 = sum(LS41),
            sumLS42 = sum(LS42),
            sumLS43 = sum(LS43),
            sumLS44 = sum(LS44),
            sumLS45 = sum(LS45),
            sumLS46 = sum(LS46),
            sumLS47 = sum(LS47),
            sumLS48 = sum(LS48),
            sumLS49 = sum(LS49),
            sumLS50 = sum(LS50),
            sumLS51 = sum(LS51),
            sumLS52 = sum(LS52),
            sumLS53 = sum(LS53),
            sumLS54 = sum(LS54),
            sumLS55 = sum(LS55),
            sumLS56 = sum(LS56),
            sumLS57 = sum(LS57),
            sumLS58 = sum(LS58),
            sumLS59 = sum(LS59),
            sumLS60 = sum(LS60),
            sumLS61 = sum(LS61),
            sumLS62 = sum(LS62),
            sumLS63 = sum(LS63),
            sumLS64 = sum(LS64),
            sumLS65 = sum(LS65),
            sumLS66 = sum(LS66),
            sumLS67 = sum(LS67),
            sumLS68 = sum(LS68),
            sumLS69 = sum(LS69),
            sumLS2 = sum(LS2),
            sumLS3 = sum(LS3),
            sumLS4 = sum(LS4),
            sumLS5 = sum(LS5),
            sumLS6 = sum(LS6),
            sumLS7 = sum(LS7),
            sumLS8 = sum(LS8),
            sumLS9 = sum(LS9))

newplantLS <- full_join(plant_data.ls, propplant.ls, by = "Facies")

newplantLS2 <- newplantLS %>%
  mutate(LS1 = LS1/sumLS1,
         LS10 = LS10/sumLS10,
         LS11 = LS11/sumLS11,
         LS12 = LS12/sumLS12,
         LS13 = LS13/sumLS13,
         LS14 = LS14/sumLS14,
         LS15 = LS15/sumLS15,
         LS16 = LS16/sumLS16,
         LS17 = LS17/sumLS17,
         LS18 = LS18/sumLS18,
         LS19 = LS19/sumLS19,
         LS20 = LS20/sumLS20,
         LS21 = LS21/sumLS21,
         LS22 = LS22/sumLS22,
         LS23 = LS23/sumLS23,
         LS24 = LS24/sumLS24,
         LS26 = LS26/sumLS26,
         LS27 = LS27/sumLS27,
         LS28 = LS28/sumLS28,
         LS29 = LS29/sumLS29,
         LS30 = LS30/sumLS30,
         LS31 = LS31/sumLS31,
         LS32 = LS32/sumLS32,
         LS33 = LS33/sumLS33,
         LS34 = LS34/sumLS34,
         LS35 = LS35/sumLS35,
         LS36 = LS36/sumLS36,
         LS37 = LS37/sumLS37,
         LS38 = LS38/sumLS38,
         LS39 = LS39/sumLS39,
         LS40 = LS40/sumLS40,
         LS41 = LS41/sumLS41,
         LS42 = LS42/sumLS42,
         LS43 = LS43/sumLS43,
         LS44 = LS44/sumLS44,
         LS45 = LS45/sumLS45,
         LS46 = LS46/sumLS46,
         LS47 = LS47/sumLS47,
         LS48 = LS48/sumLS48,
         LS49 = LS49/sumLS49,
         LS50 = LS50/sumLS50,
         LS51 = LS51/sumLS51,
         LS52 = LS52/sumLS52,
         LS53 = LS53/sumLS53,
         LS54 = LS54/sumLS54,
         LS55 = LS55/sumLS55,
         LS56 = LS56/sumLS56,
         LS57 = LS57/sumLS57,
         LS58 = LS58/sumLS58,
         LS59 = LS59/sumLS59,
         LS60 = LS60/sumLS60,
         LS61 = LS61/sumLS61,
         LS62 = LS62/sumLS62,
         LS63 = LS63/sumLS63,
         LS64 = LS64/sumLS64,
         LS65 = LS65/sumLS65,
         LS66 = LS66/sumLS66,
         LS67 = LS67/sumLS67,
         LS68 = LS68/sumLS68,
         LS69 = LS69/sumLS69,
         LS2 = LS2/sumLS2,
         LS3 = LS3/sumLS3,
         LS4 = LS4/sumLS4,
         LS5 = LS5/sumLS5,
         LS6 = LS6/sumLS6,
         LS7 = LS7/sumLS7,
         LS8 = LS8/sumLS8,
         LS9 = LS9/sumLS9)
newplantLS2[is.na(newplantLS2)] = 0

#HF/MD only plant matrix
fullplant.hfmd <- full_join(newplantHF2,newplantMD2)
fullplant.hfmd$Forest <- NULL
fullplant.hfmd$Facies <- NULL
rownames(fullplant.hfmd) <- fullplant.hfmd$Quarry
fullplant.hfmd$Quarry <- NULL
fullplant.hfmd[is.na(fullplant.hfmd)] = 0

#LS only plant matrix
#fullplant <- full_join(fullplant.hfmd,newplantLS2)
newplantLS2$Forest <- NULL
newplantLS2$Facies <- NULL
rownames(newplantLS2) <- newplantLS2$Quarry
newplantLS2$Quarry <- NULL
newplantLS2[is.na(newplantLS2)] = 0

#Plant Analysis: HF and MD----
sol.full <- metaMDS(fullplant.hfmd, distance = "bray", k=2, try = 20, trymax = 200)
sol.full
plot(sol.full, type="t") #Done with 2003 and 2014 data in dataframe
stressplot(sol.full)

#NMDS in ggplot
data.scores <- as.data.frame(scores(sol.full))  
#Using the scores function from vegan to extract the site scores 
#and convert to a data.frame
data.scores$Quarry <- rownames(data.scores) 
# create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(sol.full, "species"))  
#Using the scores function from vegan to extract the species scores 
#and convert to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
head(species.scores)  #loo

#HF and MD
data.scores$dep.env <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary",
                         "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland",
                         "Swamp", "Swamp", "Swamp", "Tributary", "Tributary",
                         "Tributary", "Fluvial", "Fluvial", "Fluvial")
data.scores$forest <- c("Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                        "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                        "SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC")

jitter <- position_jitter(width = 0.1, height = 0.1)

NMDS.plant.hfmd <- ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2)) +
  # geom_text_repel(data= species.scores,aes(label = species),
  #                 segment.size = 0.2,
  #                 point.padding = NA,
  #                 box.padding = 0.25,
  #                 nudge_y = 0.1) +
  geom_point(data=species.scores, aes(x=NMDS1,y=NMDS2), shape =4)+
  geom_point(aes(color=forest), shape = c(17,17,17,15,15,15,16,16,16,3,
                                          17,17,17,15,15,15,16,16,16), 
             size=4, position = jitter) +
  scale_color_manual(values = c("#2c6e49","#004e89"))+ 
  #guides(colour=guide_legend(reverse = TRUE))+
  #coord_equal() +
  theme_bw() +
  #ggtitle("Insect Communities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12))+
  theme(legend.position = "none")
NMDS.plant.hfmd

HFMDlegend <- get_legend(NMDS.plant.hfmd)

facies.hf <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                 "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland") 
facies.md <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                 "Tributary", "Fluvial", "Fluvial", "Fluvial") 
hfanosim <- anosim(fullplant.hfmd[c(1:10),], facies.hf, distance = "bray")
summary(hfanosim)

mdanosim <- anosim(fullplant.hfmd[c(11:19),], facies.md, distance = "bray")
summary(mdanosim)

#ANOSIM for HF and MD
facies.hfmd <- c("Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                 "Tributary", "Fluvial", "Fluvial", "Fluvial", "Upland",
                 "Swamp", "Swamp", "Swamp", "Tributary", "Tributary", 
                 "Tributary", "Fluvial", "Fluvial", "Fluvial") 

forest.hfmd <- c("Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                 "Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest","Harvard Forest",
                 "SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC","SERC")
fullplant.anosim <- anosim(fullplant.hfmd, facies.hfmd, distance = "bray") #facies a priori
summary(fullplant.anosim)

fullplant.anosim <- anosim(fullplant.hfmd, forest.hfmd, distance = "bray") #forest a piori
summary(fullplant.anosim)

#Plant Analysis: LS----
sol.full <- metaMDS(newplantLS2, distance = "bray", k=2, try = 20, trymax = 200)
sol.full
plot(sol.full, type="t") #Done with 2003 and 2014 data in dataframe
stressplot(sol.full)

#NMDS in ggplot
data.scores <- as.data.frame(scores(sol.full))  
#Using the scores function from vegan to extract the site scores 
#and convert to a data.frame
data.scores$Quarry <- rownames(data.scores) 
# create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(sol.full, "species"))  
#Using the scores function from vegan to extract the species scores 
#and convert to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores
head(species.scores)  #loo

#La Selva
data.scores$dep.env <- c("Fluvial", "Fluvial", "Fluvial", "Tributary","Tributary","Tributary", "Swamp","Swamp", "Swamp")
data.scores$forest <- c("La Selva","La Selva","La Selva","La Selva","La Selva","La Selva","La Selva","La Selva","La Selva")

NMDS.plant.ls <- ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2)) +
  # geom_text_repel(data= species.scores,aes(label = species),
  #                 segment.size = 0.2,
  #                 point.padding = NA,
  #                 box.padding = 0.25,
  #                 nudge_y = 0.1) +
  geom_point(data=species.scores, aes(x=NMDS1,y=NMDS2), shape =4)+
  geom_point(aes(color=forest), shape = c(16,16,16,15,15,15,17,17,17), 
             size=4, position = jitter) +
  scale_color_manual(values = c("#ff6b35"))+ #,"#2c6e49","#004e89"
  #guides(colour=guide_legend(reverse = TRUE))+
  #coord_equal() +
  theme_bw() +
  #ggtitle("Insect Communities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 12))+
  theme(legend.position = "none")
NMDS.plant.ls

LSlegend <- get_legend(NMDS.plant.ls)

facies.ls <- c("Fluvial","Fluvial", "Fluvial", "Tributary","Tributary","Tributary", "Swamp","Swamp", "Swamp")

fullplant.anosim <- anosim(newplantLS2, facies.ls, distance = "bray") #Can only do facies because only one forest
summary(fullplant.anosim)


#Plotting all plots together----
#all forest plots
fullfor <- plot_grid(NMDS.dt.full,
                     NMDS.plant.hfmd,
                     NMDS.plant.ls,
                     fulllegend,
                     nrow = 1,
                     ncol = 4,
                     rel_widths = c(1,1,1,0.5),
                     rel_heights = c(1,1,1,1))
fullfor

#saving----
ggsave("figures/fullforest.NMDS.pdf", fullfor)
