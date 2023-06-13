#co-occurrence analysis at the site level
library(cooccur)

#read site data from files
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

HFall_co1=(HFall_clean[,-c(1,dim(HFall_clean)[2])])
HFall_co1[HFall_co1>1]=1
HFall_co1=t(HFall_co1)


#split by quarry
HF=split(HFall_clean,HFall_clean$quarryname)
LS=split(LSall_clean,LSall_clean$quarryname)
MD=split(MDall_clean,MDall_clean$quarryname)

#combine all sites into one list
allsites=c(LS,HF,MD)

#calculate co-occurences
HF.co=cooccur(mat=HFall_co1,type = "spp_site", thresh = TRUE, true_rand_classifier = 0.1, spp_names = TRUE)

LSall_co1=LSall_clean[-c(1,dim(LSall_clean)[2])]
LSall_co1[LSall_co1>1]=1
LSall_co1=t(LSall_co1)

LS.co=cooccur(mat=LSall_co1,type = "spp_site", thresh = TRUE, true_rand_classifier = 0.1, spp_names = TRUE)

MDall_co1=MDall_clean[-c(1,dim(MDall_clean)[2])]
MDall_co1[MDall_co1>1]=1
MDall_co1=t(MDall_co1)

MD.co=cooccur(mat=MDall_co1,type = "spp_site", thresh = TRUE, true_rand_classifier = 0.1, spp_names = TRUE)



##write co-occur plot functions
#co-occurr plot functions
library(reshape2)
library(ggplot2)
paircomp=function (mod) 
{
  ptab <- mod$results
  spp <- unique(c(ptab$sp1, ptab$sp2))
  profiles <- data.frame(spp = c(), rand = c(), pos = c(), 
                         neg = c())
  for (i in as.character(spp)) {
    spptab <- ptab[ptab$sp1 == i | ptab$sp2 == i, ]
    pos <- nrow(spptab[spptab$p_gt < 0.05, ])
    neg <- nrow(spptab[spptab$p_lt < 0.05, ])
    rand <- nrow(spptab) - (pos + neg)
    total <- rand + pos + neg
    if (total == 0) {
      profiles[i, "pos"] <- 0
      profiles[i, "neg"] <- 0
      profiles[i, "rand"] <- 100
      profiles[i, "spp"] <- i
    }
    else {
      profiles[i, "pos"] <- pos/total * 100
      profiles[i, "neg"] <- neg/total * 100
      profiles[i, "rand"] <- rand/total * 100
      profiles[i, "spp"] <- i
    }
  }
  profiles$cooccur <- profiles$pos + profiles$neg
  profiles <- profiles[with(profiles, order(cooccur)), ]
  profiles$spp <- factor(x = profiles$spp, ordered = T, levels = profiles$spp)
  profiles$cooccur <- NULL
  pos <- mod$positive
  neg <- mod$negative
  rand <- mod$random
  total <- pos + neg + rand
  pos <- mod$positive/total * 100
  neg <- mod$negative/total * 100
  rand <- mod$random/total * 100
  if ("sp1_name" %in% colnames(ptab)) {
    spp_key <- unique(data.frame(sppnum = c(ptab$sp1, ptab$sp2), 
                                 sppname = c(as.character(ptab$sp1_name), as.character(ptab$sp2_name))))
    names <- merge(x = data.frame(order = 1:length(row.names(profiles)), 
                                  profiles), y = spp_key, by.x = "spp", by.y = "sppnum", 
                   all.x = T)
    names <- names[with(names, order(order)), ]
    names$order <- NULL
    names$spp <- NULL
    names$sppname <- factor(x = names$sppname, ordered = T, 
                            levels = names$sppname)
    sidebar <- data.frame(pos = pos, neg = neg, rand = rand, 
                          sppname = "All DTs")
    names <- rbind(names, sidebar)
    md <- melt(names, id = (c("sppname")))
    sppname <- md$sppname
    value <- md$value
    variable <- md$variable
    p <- ggplot(data = md, aes(x = sppname, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = "stack", 
               width = 1)
    p <- p + scale_fill_manual(values = c("greenyellow", 
                                          "lightcoral", "honeydew"), name = "", labels = c("positive", 
                                                                                           "negative", "random")) + theme(plot.title = element_text(vjust = 2, 
                                                                                                                                                    size = 20, face = "bold"), legend.text = element_text(size = 12), 
                                                                                                                          axis.title.y = element_text(size = 13), axis.text.y = element_text(size = 12), 
                                                                                                                          axis.text.x = element_text(size = 7, hjust = 0, 
                                                                                                                                                     vjust = 1, angle = -45), panel.background = element_rect(fill = "white", 
                                                                                                                                                                                                              colour = "black"), panel.grid.major = element_blank(), 
                                                                                                                          panel.grid.minor = element_blank()) + xlab("") + 
      ylab("Percent of pairings")
    p <- p + scale_y_continuous(expand = c(0.005, 
                                           0)) + scale_x_discrete(expand = c(0.005, 0))
    p <- p + geom_bar(stat = "identity", data = md[(md$sppname == 
                                                      "All DTs"), ], aes(sppname), alpha = 0, size = 0.1, 
                      color = "white", width = 1)
    p
  }
  else {
    sidebar <- data.frame(pos = pos, neg = neg, rand = rand, 
                          spp = "DTs Species")
    profiles <- rbind(profiles, sidebar)
    md <- melt(profiles, id = (c("spp")))
    p <- ggplot(data = md, aes(x = spp, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = "stack", 
               width = 1)
    p <- p + scale_fill_manual(values = c("greenyellow", 
                                          "ligtcoral", "honeydew"), name = "", labels = c("positive", 
                                                                                          "negative", "random")) + theme(plot.title = element_text(vjust = 2, 
                                                                                                                                                   size = 20, face = "bold"), legend.text = element_text(size = 18), 
                                                                                                                         axis.title.y = element_text(size = 13), axis.text.y = element_text(size = 12), 
                                                                                                                         axis.text.x = element_text(size = 12, hjust = 0, 
                                                                                                                                                    vjust = 1), panel.background = element_rect(fill = "white", 
                                                                                                                                                                                                colour = "black"), panel.grid.major = element_blank(), 
                                                                                                                         panel.grid.minor = element_blank()) + xlab("") + 
      ylab("Percent of pairings")
    p <- p + scale_y_continuous(expand = c(0.005, 
                                           0)) + scale_x_discrete(expand = c(0.005, 0))
    p <- p + geom_bar(stat = "identity", data = md[(md$spp == 
                                                      "All DTs"), ], aes(spp), alpha = 0, size = 0.1, 
                      color = "white", width = 1) + theme(axis.text.x = element_text(size=5,hjust = 0, 
                                                                                     vjust = 1, angle = -45))
    p
  }
}

#library(ggplot2)
library(ggthemes)
expobs=function (mod) 
{
  ptab <- mod$results
  ptab$signs <- ifelse(ptab$p_gt >= 0.05, 0, 1) + ifelse(ptab$p_lt >= 
                                                           0.05, 0, -1)
  exp_cooccur <- ptab$exp_cooccur
  obs_cooccur <- ptab$obs_cooccur
  signs <- ptab$signs
  p <- ggplot(ptab, aes(x = exp_cooccur, y = obs_cooccur)) + 
    geom_point(aes(fill = factor(signs, levels = c(-1, 0, 
                                                   1))),  pch = 21, size = 4)+theme_few()
  p <- p + scale_fill_manual(values = c("lightcoral", "honeydew", 
                                        "greenyellow"), name = "", labels = c("negative", "random", 
                                                                              "positive"), drop = FALSE)
  p <- p + theme(plot.title = element_text(vjust = 12, size = 15, 
                                           face = "bold"), legend.text = element_text(size = 12), 
                 axis.title = element_text(size = 15), axis.text = element_text(size = 15), 
                 axis.text.x = element_text(hjust = 0.5, vjust = 0)) + xlab("Expected Co-occurrences") + 
    ylab("Observed Co-occurrences")
  p <- p + ggtitle("Observed-Expected Plot") + geom_abline(color = "dark gray")
  p
}

coheatmap=function(cofile, dimnames){
  files=cofile[[2]]
  n=cofile$species
  vals=ifelse(files$p_gt >= 0.05, 0, 1)+ifelse(files$p_lt >= 0.05, 0, -1)
  DTmatrix=matrix(0,n,n)
  for (z in 1:nrow(files)){
    i=files$sp1[z]
    j=files$sp2[z]
    DTmatrix[i,j]=vals[z]
    DTmatrix[j,i]=DTmatrix[i,j]
  }
  
  #contdata=as.data.frame(contmatrix)
  
  colnames(DTmatrix)=dimnames
  rownames(DTmatrix)=colnames(DTmatrix)
  
  datum=as.data.frame(cbind(files$sp1,files$sp2, vals))
  
  edge_signs=datum[datum$vals!= 0,]
  
  edge_signs$col[edge_signs$vals==-1]="lightcoral"
  edge_signs$col[edge_signs$vals==1]="greenyellow"
  mat2=abs(DTmatrix)
  
  ss=which(rowSums(mat2)==0)
  
  DTmatF=DTmatrix[-ss,-ss]
  
  heatmap(DTmatF, symm = T, col = c("lightcoral", "honeydew", "yellowgreen"))
}

library(signnet)
signSBM=function(cofile, dimnames){
  files=cofile[[2]]
  n=cofile$species
  vals=ifelse(files$p_gt >= 0.05, 0, 1)+ifelse(files$p_lt >= 0.05, 0, -1)
  DTmatrix=matrix(0,n,n)
  for (z in 1:nrow(files)){
    i=files$sp1[z]
    j=files$sp2[z]
    DTmatrix[i,j]=vals[z]
    DTmatrix[j,i]=DTmatrix[i,j]
  }
  
  colnames(DTmatrix)=dimnames
  rownames(DTmatrix)=colnames(DTmatrix)
  
  datum=as.data.frame(cbind(files$sp1,files$sp2, vals))
  
  edge_signs=datum[datum$vals!= 0,]
  
  edge_signs$col[edge_signs$vals==-1]="lightcoral"
  edge_signs$col[edge_signs$vals==1]="greenyellow"
  
  g=graph_from_adjacency_matrix(abs(DTmatrix), mode = "undirected")
  
  E(g)$sign=as.numeric(edge_signs$vals)
  
  bm=signed_blockmodel(g,k = 4,alpha=0.5,annealing = TRUE)
  ggblock(g,bm$membership,show_blocks = TRUE)
  
}

paircomp2=function (mod) 
{
  ptab <- mod$results
  spp <- unique(c(ptab$sp1, ptab$sp2))
  profiles <- data.frame(spp = c(),  pos = c(), 
                         neg = c(), rand = c())
  for (i in as.character(spp)) {
    spptab <- ptab[ptab$sp1 == i | ptab$sp2 == i, ]
    pos <- nrow(spptab[spptab$p_gt < 0.05, ])
    neg <- nrow(spptab[spptab$p_lt < 0.05, ])
    rand <- nrow(spptab) - (pos + neg)
    total <- rand + pos + neg
    if (total == 0) {
      profiles[i, "rand"] <- 100
      profiles[i, "pos"] <- 0
      profiles[i, "neg"] <- 0
      profiles[i, "spp"] <- i
    }
    else {
      profiles[i, "rand"] <- rand/total * 100
      profiles[i, "pos"] <- pos/total * 100
      profiles[i, "neg"] <- neg/total * 100
      profiles[i, "spp"] <- i
    }
  }
  profiles$cooccur <- profiles$pos + profiles$neg
  profiles <- profiles[with(profiles, order(cooccur)), ]
  profiles$spp <- factor(x = profiles$spp, ordered = T, levels = profiles$spp)
  profiles$cooccur <- NULL
  pos <- mod$positive
  neg <- mod$negative
  rand <- mod$random
  total <- pos + neg + rand
  pos <- mod$positive/total * 100
  neg <- mod$negative/total * 100
  rand <- mod$random/total * 100
  if ("sp1_name" %in% colnames(ptab)) {
    spp_key <- unique(data.frame(sppnum = c(ptab$sp1, ptab$sp2), 
                                 sppname = c(as.character(ptab$sp1_name), as.character(ptab$sp2_name))))
    names <- merge(x = data.frame(order = 1:length(row.names(profiles)), 
                                  profiles), y = spp_key, by.x = "spp", by.y = "sppnum", 
                   all.x = T)
    names <- names[with(names, order(order)), ]
    names$order <- NULL
    names$spp <- NULL
    names$sppname <- factor(x = names$sppname, ordered = T, 
                            levels = names$sppname)
    sidebar <- data.frame(rand = rand, pos = pos, neg = neg, 
                          sppname = "All DTs")
    names <- rbind(names, sidebar)
    md <- melt(names, id = (c("sppname")))
    sppname <- md$sppname
    value <- md$value
    variable <- md$variable
    p <- ggplot(data = md, aes(x = sppname, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = "stack", 
               width = 1)
    p <- p + scale_fill_manual(values = c("honeydew","greenyellow", 
                                          "lightcoral"), name = "", labels = c("random","positive", 
                                                                               "negative")) + theme(plot.title = element_text(vjust = 2, 
                                                                                                                              size = 20, face = "bold"), legend.text = element_text(size = 12), 
                                                                                                    axis.title.y = element_text(size = 13), axis.text.y = element_text(size = 12), 
                                                                                                    axis.text.x = element_text(size = 7, hjust = 0, 
                                                                                                                               vjust = 1, angle = -45), panel.background = element_rect(fill = "white", 
                                                                                                                                                                                        colour = "black"), panel.grid.major = element_blank(), 
                                                                                                    panel.grid.minor = element_blank()) + xlab("") + 
      ylab("Percent of pairings")
    p <- p + scale_y_continuous(expand = c(0.005, 
                                           0)) + scale_x_discrete(expand = c(0.005, 0))
    p <- p + geom_bar(stat = "identity", data = md[(md$sppname == 
                                                      "All DTs"), ], aes(sppname), alpha = 0, size = 0.1, 
                      color = "white", width = 1)
    p
  }
  else {
    sidebar <- data.frame(pos = pos, neg = neg, rand = rand, 
                          spp = "DTs Species")
    profiles <- rbind(profiles, sidebar)
    md <- melt(profiles, id = (c("spp")))
    p <- ggplot(data = md, aes(x = spp, y = value, fill = variable)) + 
      geom_bar(stat = "identity", position = "stack", 
               width = 1)
    p <- p + scale_fill_manual(values = c( "honeydew","greenyellow", 
                                           "ligtcoral"), name = "", labels = c("random","positive", 
                                                                               "negative")) + theme(plot.title = element_text(vjust = 2, 
                                                                                                                              size = 20, face = "bold"), legend.text = element_text(size = 18), 
                                                                                                    axis.title.y = element_text(size = 13), axis.text.y = element_text(size = 12), 
                                                                                                    axis.text.x = element_text(size = 12, hjust = 0, 
                                                                                                                               vjust = 1), panel.background = element_rect(fill = "white", 
                                                                                                                                                                           colour = "black"), panel.grid.major = element_blank(), 
                                                                                                    panel.grid.minor = element_blank()) + xlab("") + 
      ylab("Percent of pairings")
    p <- p + scale_y_continuous(expand = c(0.005, 
                                           0)) + scale_x_discrete(expand = c(0.005, 0))
    p <- p + geom_bar(stat = "identity", data = md[(md$spp == 
                                                      "All DTs"), ], aes(spp), alpha = 0, size = 0.1, 
                      color = "white", width = 1) + theme(axis.text.x = element_text(size=5,hjust = 0, 
                                                                                     vjust = 1, angle = -45))
    p
  }
}

#plot them
expobs(HF.co)

#export 8x6
paircomp(HF.co)
paircomp(LS.co)
paircomp(MD.co)

paircomp2(HF.co)
paircomp2(LS.co)
paircomp2(MD.co)

#export as 8x8 cm
coheatmap(HF.co, rownames(HFall_co1))
coheatmap(LS.co, rownames(LSall_co1))
coheatmap(MD.co, rownames(MDall_co1))
