
#read three site data
HFall=read.csv("data/abundance_data_HFmatrix2.csv")
LSall=read.csv("data/abundance_data_LSmatrix2.csv")
MDall=read.csv("data/abundance_data_SERCmatrix2.csv")

#quarry names
HFnames=paste("HF", unique(HFall$quarry))
LSnames=paste("LS", unique(HFall$quarry))
MDnames=paste("HF", unique(HFall$quarry))

HFall$quarryname=paste(HFall$site,".", HFall$quarry, sep="")
LSall$quarryname=paste(LSall$site,".", LSall$quarry, sep="")
MDall$quarryname=paste(MDall$site,".", MDall$quarry, sep="")

#clean the other unreq columns and remove None DT
HFall_clean=HFall[,-c(1:5,7:17,26)]
colnames(HFall_clean)[233:235]=c("DT.223", "DT.261", "DT.276")
LSall_clean=LSall[,-c(1:5,7:17,25)]
MDall_clean=MDall[,-c(1:5,7:17,27)]

#split by quarry
HF=split(HFall_clean,HFall_clean$quarryname)
LS=split(LSall_clean,LSall_clean$quarryname)
MD=split(MDall_clean,MDall_clean$quarryname)

#remove the sites with less than 300 leaves
LS=LS[-c(5,6)]

#combine all sites into one list
allsites=c(LS,HF,MD)

set.seed(1)

networklevlist=NULL
networklevmean=NULL
networklevelvar=NULL

library(bipartite)
library(Rfast)
library(dplyr)

for (i in 1:length(allsites)){
  nam=names(allsites)[[i]]
  netcurr=allsites[[i]]
  netcurr2=netcurr %>% group_by(species) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
  netmat = data.matrix(netcurr2[,-c(1)], rownames.force = NA)
  rownames(netmat)= netcurr2$species
  
  minrow=300
  
  netpro=NULL
  
  allnames=cbind(rownames(netmat),0)
  colnames(allnames)=c("species","Nul")
  allnames=as.data.frame(allnames)
  
  for (j in 1:100){
    q2=sample(nrow(netcurr), minrow, replace=F)
    datarand=netcurr[q2,]
    netrand=datarand %>% group_by(species) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
    netrandfinal=left_join(allnames,netrand, by="species")
    netrandfinal[is.na(netrandfinal)] <- 0
    netrandfinal=netrandfinal[,-c(1:2)]
    netrandfinal=data.matrix(netrandfinal)
    netpro=rbind(netpro,networklevel(netrandfinal))
    
    if (j %% 20 == 0){print(paste("Sequence", " ", j))}
  }
  
  
  
  nam2=nam
  
  networklevmean=rbind(networklevmean,c(nam2,colmeans(netpro)))
  networklevelvar=rbind(networklevelvar,c(nam2,colVars(netpro, std = F)))
  
  networklevlist[[i]]=netpro
  
  print(i)
}

colnames(networklevmean)=c("Web",colnames(netpro))
colnames(networklevelvar)=c("Web",colnames(netpro))

networklevmean=read.csv("results/bootstrapmean.csv")


library(dplyr)
library(bipartite)
networklev=NULL
higherlev=NULL
plantlev=NULL
plantnam=NULL
highernam=NULL

for (i in 1:length(allsites)){
  nam=names(allsites)[[i]]
  netcurr=allsites[[i]]
  netcurr2=netcurr %>% group_by(species) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
  netmat = data.matrix(netcurr2[,-c(1)], rownames.force = NA)
  rownames(netmat)= netcurr2$species
  nam2=nam
  #plotweb(net)
  networklev=rbind(networklev,c(nam2,networklevel(netmat)))
  X=specieslevel(netmat)
  X1=X$`higher level`
  X2=X$`lower level`
  
  X1$web=nam2
  X2$web=nam2
  plantnam=c(plantnam,rownames(X2))
  highernam=c(highernam,rownames(X1))
  
  higherlev=rbind(higherlev,X1)
  plantlev=rbind(plantlev,X2)
  print(i)
}


plantlev$speciesID=plantnam
higherlev$speciesID=highernam

higherlev$site="NN"
for (i in 1:nrow(higherlev)) {
  higherlev$site[i]=strsplit(as.character(higherlev$web), split =  "\\.")[[i]][1]
}

plantlev$site="NN"
for (i in 1:nrow(plantlev)) {
  plantlev$site[i]=strsplit(as.character(plantlev$web), split =  "\\.")[[i]][1]
}
