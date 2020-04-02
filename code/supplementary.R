#######################################################################################
#
#                                                                                     
#   Filename    :	      supplementary.R (UTF-8)  												  
#                                                                                     
#   Project     :       bioRxiv preprint "Tracking the outbreak. 
#                       An optimized delimiting survey strategy for Xylella fastidiosa
#                       in Alicante, Spain"                                                             
#   Authors     :       E. Lázaro et al.                                                              
#   Date        :       01.03.2020
#   Purpose     :       Reproduce supplementary material 
#																				  
#                                                 
########################################################################################






# Packages
library(foreign)
library(INLA)
library(rgdal)


# Table S1

# Load selected models. WAIC and LCPO selected models

load("./results/incidence/models/model_selection.RData")

waic <- NULL
lcpo <- NULL


for(i in 1:16){
  
  model = load(paste('./results/incidence/models/mod',i,'.RData', sep=""))
  
  models_lcpo = load(paste('./results/incidence/models/mod',i,'_lcpo','.RData', sep=""))
  
  waic[i] = round(get(paste("mod",i,sep=""))$waic$waic, dig=3)
  
  lcpo[i] = round(-mean(log(get(paste("mod",i,"_lcpo", sep=""))$cpo$cpo)), dig=3)
}

tables1 <- cbind(as.character(model_selection$`Models according to waic`$Models[1:16]),as.numeric(waic),as.numeric(lcpo))
colnames(tables1) <- c("Model","WAIC","LCPO")
tables1

###########################################################################

# Table S2

load("./data/binomial_subsets/data1_9.RData") #DS9
load("./data/binomial_subsets/data1_23.RData") #DS23
load("./data/binomial_subsets/data1_37.RData") #DS37
load("./data/binomial_subsets/data1_51.RData") #DS51
load("./data/binomial/data1.RData")# whole database


replicates = 100
nsamples_9 <- list()
nsamples_23 <- list()
nsamples_37 <-list()
nsamples_51 <-list()

npositives_9 <- list()
npositives_23 <- list()
npositives_37 <-list()
npositives_51 <-list()

for(i in 1:replicates){
  
  nsamples_9[[i]]<-as.vector(sum(data1_9[[i]]$nsamples))
  nsamples_23[[i]]<-as.vector(sum(data1_23[[i]]$nsamples))
  nsamples_37[[i]]<-as.vector(sum(data1_37[[i]]$nsamples))
  nsamples_51[[i]]<-as.vector(sum(data1_51[[i]]$nsamples))
  
  npositives_9[[i]]<-as.vector(sum(data1_9[[i]]$positives))
  npositives_23[[i]]<-as.vector(sum(data1_23[[i]]$positives))
  npositives_37[[i]]<-as.vector(sum(data1_37[[i]]$positives))
  npositives_51[[i]]<-as.vector(sum(data1_51[[i]]$positives))
}


ds9 <- c(median(unlist(nsamples_9)), median(unlist(npositives_9)), -median(unlist(npositives_9))+ median(unlist(nsamples_9)),median(unlist(npositives_9))/median(unlist(nsamples_9)))
ds9 <- round(ds9, dig=3)

ds23 <- c(median(unlist(nsamples_23)), median(unlist(npositives_23)), -median(unlist(npositives_23))+ median(unlist(nsamples_23)),median(unlist(npositives_23))/median(unlist(nsamples_23)))
ds23 <- round(ds23, dig=3)

ds37 <- c(median(unlist(nsamples_37)), median(unlist(npositives_37)), -median(unlist(npositives_37))+ median(unlist(nsamples_37)),median(unlist(npositives_37))/median(unlist(nsamples_37)))
ds37 <- round(ds37, dig=3)

ds51 <- c(median(unlist(nsamples_51)), median(unlist(npositives_51)), -median(unlist(npositives_51))+ median(unlist(nsamples_51)),median(unlist(npositives_51))/median(unlist(nsamples_51)))
ds51 <- round(ds51, dig=3)

reference <- c(sum(data1$nsamples), sum(data1$positives), -sum(data1$positives)+sum(data1$nsamples), sum(data1$positives)/sum(data1$nsamples) )
reference <- round(reference, dig=3)


tables2 <- rbind(ds9,ds23,ds37,ds51,reference)
colnames(tables2) <- c("Nsamples", "Npositives","Nnegatives","Global incidence")
tables2


# Figures S1 and S2

Grid1_2018 <- readOGR("./data/grids", "Grid1_2018") # 1km x 1km grid (whole survey area)


# DS9
# Sampling intenisty

nsamples <- data1_9[[which(unlist(npositives_9)==median(unlist(npositives_9)))[1]]]$nsamples

nsamplesf <- cut(nsamples, 
                 breaks = c(1,2,5,9,23,37,51, 109, Inf ), 
                 labels = c("1","[2,4]","[5,8]", "[9,22]", "[23,36]", "[37,50]","[51,108]",+109), 
                 right = FALSE)

Grid1_2018$nsamplesf <- nsamplesf

# N positives/cell

positives <- data1_9[[which(unlist(npositives_9)==median(unlist(npositives_9)))[1]]]$positives

positivesf <- cut(positives, 
                  breaks = c(-Inf, 1, 2, 3, 5,8,12, Inf ), 
                  labels = c(0, "1", "2 ","3", "[5,7]","[8,11]","+12"), 
                  right = FALSE)

Grid1_2018$positivesf <- positivesf

# Presence absence cells

Grid1_2018$Res <- as.factor(ifelse(data1_9[[which(unlist(npositives_9)==median(unlist(npositives_9)))[1]]]$positives>=1,1,0))

levels(Grid1_2018$Res) <- c("Absence","Presence")



# Incicende per cell

incidence <- positives/nsamples

incidencef <- cut(incidence, 
                     breaks = c(-Inf,0,0.01, 0.02, 0.05, 0.1, 0.6, +Inf ), 
                     labels = c(0, "]0,0.01]","]0.01,0.02] ", "]0.02,0.05] ","]0.05,0.1]", "]0.1,0.6]","]0.6,1]"), 
                     right = FALSE)

Grid1_2018$incidencef <- incidencef




paleta <- colorRampPalette(c("Yellow2","Orange2","Brown4"))(10)

jpeg(filename="./s1a1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0),oma=c(0.001,0.001,0.001,0.001) ,cex.axis=10, cex.lab=10, cex.main=10)
spplot(Grid1_2018, "nsamplesf",cex=10,col.regions =paleta, colorkey=F)

#box(col = 'black')
dev.off()

jpeg(filename="./s1b1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("positivesf"), names.attr = c("Total"),
       col.regions =paleta,
       colorkey=F)
#box(col = 'black')
dev.off()

jpeg(filename="./s2a1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Res"), names.attr = c("Total"),col.regions=c("Yellow2","Brown4"),
       colorkey=F)
#box(col = 'black')
dev.off()




jpeg(filename="./s2b1.jpeg",  
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("incidencef"), names.attr = c("Total"),col.regions =paleta,
       colorkey=F)
dev.off()

# DS23
# Sampling intenisty

nsamples <- data1_23[[which(unlist(npositives_23)==median(unlist(npositives_23)))[1]]]$nsamples

nsamplesf <- cut(nsamples, 
                 breaks = c(1,2,5,9,23,37,51, 109, Inf ), 
                 labels = c("1","[2,4]","[5,8]", "[9,22]", "[23,36]", "[37,50]","[51,108]",+109), 
                 right = FALSE)

Grid1_2018$nsamplesf <- nsamplesf

# N positives/cell

positives <- data1_23[[which(unlist(npositives_23)==median(unlist(npositives_23)))[1]]]$positives

positivesf <- cut(positives, 
                  breaks = c(-Inf, 1, 2, 3, 5,8,12, Inf ), 
                  labels = c(0, "1", "2 ","3", "[5,7]","[8,11]","+12"), 
                  right = FALSE)

Grid1_2018$positivesf <- positivesf

# Presence absence cells

Grid1_2018$Res <- as.factor(ifelse(data1_23[[which(unlist(npositives_23)==median(unlist(npositives_23)))[1]]]$positives>=1,1,0))

levels(Grid1_2018$Res) <- c("Absence","Presence")



# Incicende per cell

incidence <- positives/nsamples

incidencef <- cut(incidence, 
                  breaks = c(-Inf,0,0.01, 0.02, 0.05, 0.1, 0.6, +Inf ), 
                  labels = c(0, "]0,0.01]","]0.01,0.02] ", "]0.02,0.05] ","]0.05,0.1]", "]0.1,0.6]","]0.6,1]"), 
                  right = FALSE)

Grid1_2018$incidencef <- incidencef




paleta <- colorRampPalette(c("Yellow2","Orange2","Brown4"))(10)

jpeg(filename="./s1a2.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0.001,0.001,0.001,0.001) ,cex.axis=10, cex.lab=10, cex.main=10)
spplot(Grid1_2018, "nsamplesf",cex=10,col.regions =paleta, colorkey=F)

#box(col = 'black')
dev.off()

jpeg(filename="./s1b2.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("positivesf"), names.attr = c("Total"),
       col.regions =paleta,
       colorkey=F)
#box(col = 'black')
dev.off()

jpeg(filename="./s2a2.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Res"), names.attr = c("Total"),col.regions=c("Yellow2","Brown4"),
       colorkey=F)
#box(col = 'black')
dev.off()




jpeg(filename="./s2b2.jpeg",  
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("incidencef"), names.attr = c("Total"),col.regions =paleta,
       colorkey=F)
dev.off()



# DS37
# Sampling intenisty

nsamples <- data1_37[[which(unlist(npositives_37)==median(unlist(npositives_37)))[1]]]$nsamples

nsamplesf <- cut(nsamples, 
                 breaks = c(1,2,5,9,23,37,51, 109, Inf ), 
                 labels = c("1","[2,4]","[5,8]", "[9,22]", "[23,36]", "[37,50]","[51,108]",+109), 
                 right = FALSE)

Grid1_2018$nsamplesf <- nsamplesf

# N positives/cell

positives <- data1_37[[which(unlist(npositives_37)==median(unlist(npositives_37)))[1]]]$positives

positivesf <- cut(positives, 
                  breaks = c(-Inf, 1, 2, 3, 5,8,12, Inf ), 
                  labels = c(0, "1", "2 ","3", "[5,7]","[8,11]","+12"), 
                  right = FALSE)

Grid1_2018$positivesf <- positivesf

# Presence absence cells

Grid1_2018$Res <- as.factor(ifelse(data1_37[[which(unlist(npositives_37)==median(unlist(npositives_37)))[1]]]$positives>=1,1,0))

levels(Grid1_2018$Res) <- c("Absence","Presence")



# Incicende per cell

incidence <- positives/nsamples

incidencef <- cut(incidence, 
                  breaks = c(-Inf,0,0.01, 0.02, 0.05, 0.1, 0.6, +Inf ), 
                  labels = c(0, "]0,0.01]","]0.01,0.02] ", "]0.02,0.05] ","]0.05,0.1]", "]0.1,0.6]","]0.6,1]"), 
                  right = FALSE)

Grid1_2018$incidencef <- incidencef




paleta <- colorRampPalette(c("Yellow2","Orange2","Brown4"))(10)

jpeg(filename="./s1a3.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0.001,0.001,0.001,0.001) ,cex.axis=10, cex.lab=10, cex.main=10)
spplot(Grid1_2018, "nsamplesf",cex=10,col.regions =paleta, colorkey=F)

#box(col = 'black')
dev.off()

jpeg(filename="./s1b3.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("positivesf"), names.attr = c("Total"),
       col.regions =paleta,
       colorkey=F)
#box(col = 'black')
dev.off()

jpeg(filename="./s2a3.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Res"), names.attr = c("Total"),col.regions=c("Yellow2","Brown4"),
       colorkey=F)
#box(col = 'black')
dev.off()

jpeg(filename="./s2b3.jpeg",  
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("incidencef"), names.attr = c("Total"),col.regions =paleta,
       colorkey=F)
dev.off()

# DS51
# Sampling intenisty

nsamples <- data1_51[[which(unlist(npositives_51)==median(unlist(npositives_51)))[1]]]$nsamples

nsamplesf <- cut(nsamples, 
                 breaks = c(1,2,5,9,23,37,51, 109, Inf ), 
                 labels = c("1","[2,4]","[5,8]", "[9,22]", "[23,36]", "[37,50]","[51,108]",+109), 
                 right = FALSE)

Grid1_2018$nsamplesf <- nsamplesf

# N positives/cell

positives <- data1_51[[which(unlist(npositives_51)==median(unlist(npositives_51)))[1]]]$positives

positivesf <- cut(positives, 
                  breaks = c(-Inf, 1, 2, 3, 5,8,12, Inf ), 
                  labels = c(0, "1", "2 ","3", "[5,7]","[8,11]","+12"), 
                  right = FALSE)

Grid1_2018$positivesf <- positivesf

# Presence absence cells

Grid1_2018$Res <- as.factor(ifelse(data1_51[[which(unlist(npositives_51)==median(unlist(npositives_51)))[1]]]$positives>=1,1,0))

levels(Grid1_2018$Res) <- c("Absence","Presence")



# Incicende per cell

incidence <- positives/nsamples

incidencef <- cut(incidence, 
                  breaks = c(-Inf,0,0.01, 0.02, 0.05, 0.1, 0.6, +Inf ), 
                  labels = c(0, "]0,0.01]","]0.01,0.02] ", "]0.02,0.05] ","]0.05,0.1]", "]0.1,0.6]","]0.6,1]"), 
                  right = FALSE)

Grid1_2018$incidencef <- incidencef




paleta <- colorRampPalette(c("Yellow2","Orange2","Brown4"))(10)

jpeg(filename="./s1a4.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0.001,0.001,0.001,0.001) ,cex.axis=10, cex.lab=10, cex.main=10)
spplot(Grid1_2018, "nsamplesf",cex=10,col.regions =paleta,  colorkey=list(space="bottom", height = 1),par.settings=list(fontsize=list(text=20)))

#box(col = 'black')
dev.off()

jpeg(filename="./s1b4.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("positivesf"), names.attr = c("Total"),
       col.regions =paleta,
       colorkey=list(space="bottom", height = 1),par.settings=list(fontsize=list(text=20)))
#box(col = 'black')
dev.off()

jpeg(filename="./s2a4.jpeg", 
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Res"), names.attr = c("Total"),col.regions=c("Yellow2","Brown4"),
       colorkey=list(space="bottom", height = 1),par.settings=list(fontsize=list(text=20)))
#box(col = 'black')
dev.off()




jpeg(filename="./s2b4.jpeg",  
     units="in", 
     width=8, 
     height=6.5, 
     pointsize=12, 
     res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("incidencef"), names.attr = c("Total"),col.regions =paleta,
       colorkey=list(space="bottom", height = 1),par.settings=list(fontsize=list(text=20)))
dev.off()

