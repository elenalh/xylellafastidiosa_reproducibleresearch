#######################################################################################
#
#                                                                                     
#   Filename    :	      article.R  (UTF-8)  												  
#                                                                                     
#   Project     :       bioRxiv preprint "Tracking the outbreak. 
#                       An optimized delimiting survey strategy for Xylella fastidiosa
#                       in Alicante, Spain"                                                             
#   Authors     :       E. Lázaro et al.                                                              
#   Date        :       01.03.2020
#   Purpose     :       Reproduce main document results as described in bioRxiv preprint
#																				  
#                                                 
########################################################################################

# Working directory

setwd("~./xylellafastidiosa_reproducibleresearch")


# Packages

library(raster)
library(rgdal)
library(rgeos)
library(prettymapr)
library(RColorBrewer)
library(INLA)





# Data

# Survey data

load("./data/survey_data/xf2018.RData") # Data of 2018 delimiting Xf survey campaign in Alicante region


# Survey area boundaries

zd_2018 <- readOGR("./data/boundaries", "zd_2018") #survey area boundaries

infzone_2018 <- readOGR("./data/boundaries", "infzone_2018") # Infected zone boundaries

zone100_2018 <- readOGR("./data/boundaries", "zone100_2018") # 100m x 100m grid boundaries

# Survey grids


Grid100_2018 <- readOGR("./data/grids", "Grid100_2018") # 100m x 100 m grid (whole survey area)

Grid100_2018b <- readOGR("./data/grids", "Grid100_2018b") # 100m x 100m grid (first kilometer)

Grid1_2018 <- readOGR("./data/grids", "Grid1_2018") # 1km x 1km grid (whole survey area)

Grid500_2018 <- readOGR("./data/grids", "Grid500_2018") # 500m x 500m grid (whole survey area)


###########################################################


# 2. Material and methods 

# 2.1 Database

# Table 1

tab1 <- addmargins(table(xf2018@data$species,xf2018@data$reslab,xf2018@data$symptoms), 2)
tab1[,c(1,3),1]
tab1[,c(1,3),2]

tab2 <- addmargins(table(xf2018@data$species,xf2018@data$reslab), 2)
tab2[,c(1,3)]

table1 <- cbind(tab1[,c(1,3),1],tab1[,c(1,3),2],tab2[,c(1,3)])

#######################################################################################

# 2.2 Evaluation of delimiting strategies

# Current delimiting strategy

# Figure 1a

jpeg("./1a.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
layout(matrix(c(1, 1,2, 2), nrow=2, byrow=TRUE),heights = c(2.2,0.35))
par(mar=c(0,0.8,0,0))
plot(zd_2018)
plot(Grid1_2018,border="blue",add=T)
plot(Grid100_2018,border="darkorange",add=T)
plot(Grid100_2018b,border="darkorange",add=T)
plot(infzone_2018,col="red",add=T)
plot(Grid1_2018,border="blue",add=T)
axis(side=1,at=seq(727000,780000,by=3000), pos=4270700)
axis(side=2,at=seq(4271000,4309000,by=3000), pos= 726500)
addscalebar(pos = "bottomright", plotepsg= 32630, padin=c(0.6, 1))
par(mar=c(0.8,0,0,0))
plot(crop(Grid100_2018b,Grid1_2018[225,]),border="darkorange")
plot(Grid1_2018[225,],border="blue",add=T)
leg <- c(expression(paste("1 km"^"2")),expression(paste("0.01 km"^"2")))
legend("right" ,leg , col=c("blue","darkorange"),
       pch=c(0,0,0),bty="n",cex=1.65)
dev.off()

jpeg("./1b.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
xf2018p <- xf2018[xf2018$reslab=="positive",]
layout(matrix(c(1, 1,2, 2), nrow=2, byrow=TRUE),heights = c(2.2,0.35))
par(mar=c(0,0.8,0,0))
plot(zd_2018)
points(xf2018[xf2018$reslab=="negative",],col="green",pch=16)
points(xf2018p,col="red",pch=16)
plot(Grid1_2018,border="blue",add=T)
axis(side=1,at=seq(727000,780000,by=3000), pos=4270700)
addscalebar(pos = "bottomright", plotepsg= 32630, padin=c(0.6, 1))
par(mar=c(0.8,0,0,0))
plot(crop(Grid100_2018b,Grid1_2018[225,]),border="white")
plot(Grid1_2018[225,],border="white",add=T)
leg <- c("Positive","Negative")
legend("right" ,leg , col=c("red","green"),
       pch=c(16,16),bty="n",cex=1.65)
dev.off()

# Table 2

load("./data/binomial/data1.RData") # Binomial data grid 1 x 1
load("./data/binomial/data500.RData") # Binomial data grid 500 x 500
load("./data/binomial/data100.RData") # Binomial data grid 100 x 100

grid1 <- cbind(nrow(data1), sum(ifelse(data1$nsamples >= 1 , 1, 0)),
               sum(ifelse(data1$positives >= 1 , 1, 0)), median(data1$nsamples), max(data1$nsamples))
grid500 <- cbind(nrow(data500), sum(ifelse(data500$nsamples[!is.na(data500$positives)] >= 1 , 1, 0)),
                 sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0)),
                 median(data500$nsamples[!is.na(data500$positives)]),
                 max(data500$nsamples[!is.na(data500$positives)]))


grid100 <- cbind(nrow(data100), sum(ifelse(data100$nsamples[!is.na(data100$positives)] >= 1 , 1, 0)),
                 sum(ifelse((data100$positives >= 1) &  (!is.na(data100$positives)) * 1, 1, 0)),
                 median(data100$nsamples[!is.na(data100$positives)]),
                 max(data100$nsamples[!is.na(data100$positives)]))

table2 <- matrix(c(grid1, grid500, grid100), nrow=3,ncol=5, byrow=TRUE)
rownames(table2) <- c("j=1","j=500", "j=100")
colnames(table2) <- c("Cj","Cj,s", "Cj,+","median", "max.")
table2

# Sequential adaptive strategy

# Figure 2

# Figure 2a (3 phase design-phase1)

grida <- Grid1_2018[c(367),]
gridb <- Grid1_2018[c(368),]
set.seed(125)
pointsa <- spsample(grida,n=50,"random")
pointsb <- spsample(gridb,n=50,"random")
set.seed(125)
resulta <- sample(c(0,1),50,prob = c(0.95,0.05),replace=T)
resultb <- rep(0,50)
pointsa$result <- resulta
pointsb$result <- resultb
grid <- rbind(grida,gridb)

jpeg("./2a.jpeg", units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0))
plot(grid[,],border="red",lwd=4)
plot(grid[2,],border="green",lwd=4,add=TRUE)
plot(grid[1,],border="red",lwd=4,add=TRUE)
points(pointsa[pointsa$result=="1",], col="red", pch=16,lwd=4)
points(pointsa[pointsa$result=="0",], col="green", pch=16, lwd=4)
points(pointsb[pointsb$result=="0",], col="green", pch=16,lwd=4)
text(x=751200,y=4287100,"n1_1",cex=2)
text(x=752200,y=4287100,"n1_2",cex=2)
dev.off()


# Figure 2b  (3 phase design-phase2)

grid500_1 <- crop(Grid500_2018,grid[1,])
datos_1 <- crop(pointsa,extent<-extent(grid500_1))

jpeg("./2b.jpeg",units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0))
plot(grid,border="gray",lwd=4)
plot(grid500_1,border="green",lwd=4,add=T)
plot(grid500_1[c(2,4),],border="red",lwd=4,add=T)
points(datos_1[datos_1$result=="1",], col="red", pch=16,lwd=4)
points(datos_1[datos_1$result=="0",], col="green", pch=16,lwd=4)
text(x=751200,y=4287550,"n0.25_1",cex=2)
text(x=751700,y=4287550,"n0.25_2",cex=2)
text(x=751200,y=4287050,"n0.25_3",cex=2)
text(x=751700,y=4287050,"n0.25_4",cex=2)
dev.off()

# Figure 2c (3 phase design-phase3)

grid100_1 <- crop(Grid100_2018,grid500_1[c(2,4),])
datos_2 <- crop(datos_1,extent<-extent(grid100_1))

jpeg("./2c.jpeg",units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0))
plot(grid,border="gray",lwd=4)
plot(grid500_1,add=T,lwd=4,border="gray")
plot(grid100_1,border="green",lwd=4,add=T)
plot(grid100_1[c(15,35,40),],border="red",lwd=3,add=T)
points(datos_2[datos_2$result=="1",], col="red", pch=16)
points(datos_2[datos_2$result=="0",], col="green", pch=16)

text(x=extent(grid100_1[1,])[1],y=extent(grid100_1[1,])[3]+50,"n0.01_1",cex=2)
text(x=extent(grid100_1[50,])[2]+50,y=extent(grid100_1[50,])[3]+50,"n0.01_50",cex=2)
dev.off()


# Figure 2d (2 phase design-phase1)

grida <- Grid1_2018[c(367),]
gridb <- Grid1_2018[c(368),]
set.seed(125)
pointsa <- spsample(grida,n=50,"random")
pointsb <- spsample(gridb,n=50,"random")
set.seed(125)
resulta <- sample(c(0,1),50,prob = c(0.95,0.05),replace=T)
resultb <- rep(0,50)
pointsa$result <- resulta
pointsb$result <- resultb

grid<-rbind(grida,gridb)

jpeg("./2d.jpeg", units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0))
plot(grid[,],border="red",lwd=4)
plot(grid[2,],border="green",lwd=4,add=TRUE)
plot(grid[1,],border="red",lwd=4,add=TRUE)
points(pointsa[pointsa$result=="1",], col="red", pch=16,lwd=4)
points(pointsa[pointsa$result=="0",], col="green", pch=16, lwd=4)
points(pointsb[pointsb$result=="0",], col="green", pch=16,lwd=4)
text(x=751200,y=4287100,"n1_1",cex=2)
text(x=752200,y=4287100,"n1_2",cex=2)
dev.off()




# Figure 2f (2 phase design-phase2)

grid100_2<-crop(Grid100_2018,grid[1,])
datos_3<-crop(pointsa,extent<-extent(grid100_2))

jpeg("./2f.jpeg",units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0))
plot(grid,border="gray",lwd=4)
plot(grid100_2,border="green",lwd=4,add=T)
plot(grid100_2[c(30,70,80),],border="red",lwd=4,add=T)
points(datos_3[datos_3$result=="1",], col="red", pch=16,lwd=4)
points(datos_3[datos_3$result=="0",], col="green", pch=16,lwd=4)
text(x=extent(grid100_2[1,])[1]+160,y=extent(grid100_2[1,])[3]+50,"n0.01_1",cex=2)
text(x=extent(grid100_2[10,])[2]+160,y=extent(grid100_2[10,])[3]+50,"n0.01_10",cex=2)
text(x=extent(grid100_2[91,])[1]+160,y=extent(grid100_2[91,])[3]+50,"n0.01_91",cex=2)
text(x=extent(grid100_2[100,])[2]+160,y=extent(grid100_2[100,])[3]+50,"n0.01_100",cex=2)
dev.off()



########################################################

# 3. Results

# 3.1 Evaluation of delimiting strategies

# Table 3

jpeg("./t3a.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
xf2018p <- xf2018[xf2018$ResLab=="POSITIVO",]
layout(matrix(c(1, 1,2, 2), nrow=2, byrow=TRUE),heights = c(2.2,0.35))
par(mar=c(0,0.8,0,0))
plot(zd_2018)
plot(Grid1_2018,border="blue",add=T)
Grid1_2018 <- spTransform(Grid1_2018,CRSobj = crs(xf2018p))
aux <- over(Grid1_2018,xf2018p)
aux <- crop(Grid500_2018,Grid1_2018[-which(is.na(aux$X)),])
aux <- spTransform(aux,CRSobj = crs(xf2018p))
plot(aux, border="cornsilk4",add=T)
axis(side=2,at=seq(4271000,4309000,by=3000), pos= 726500)
addscalebar(pos = "bottomright", plotepsg= 32630, padin=c(0.6, 1))
Grid100_2018b <- spTransform(Grid100_2018b,CRSobj = crs(xf2018p))
aux2 <- over(aux,xf2018p)
plot(crop(Grid100_2018,aux[-which(is.na(aux2$X)),]), border="darkorange", add=T)
plot(Grid1_2018,border="blue",add=T)
par(mar=c(0.8,0,0,0))
plot(crop(Grid100_2018,Grid1_2018[225,]),border="darkorange")
plot(crop(Grid500_2018,Grid1_2018[225,]),border="cornsilk4", add=T)
plot(Grid1_2018[225,],border="blue",add=T)
leg <- c(expression(paste("1 km"^"2")), expression(paste("0.25 km"^"2")), expression(paste("0.01 km"^"2")))
legend("right" ,leg , col=c("blue","cornsilk4","darkorange"),
       pch=c(0,0,0),bty="n",cex=1.65)
dev.off()

jpeg("./t3b.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
layout(matrix(c(1, 1,2, 2), nrow=2, byrow=TRUE),heights = c(2.2,0.35))
par(mar=c(0,0.8,0,0))
plot(zd_2018)
plot(Grid1_2018,border="blue",add=T)
Grid1_2018 <- spTransform(Grid1_2018,CRSobj = crs(xf2018p))
aux <- over(Grid1_2018,xf2018p)
plot(crop(Grid100_2018,Grid1_2018[-which(is.na(aux$X)),]),border="darkorange",add=T)
plot(Grid1_2018,border="blue",add=T)
axis(side=2,at=seq(4271000,4309000,by=3000), pos= 726500)
addscalebar(pos = "bottomright", plotepsg= 32630, padin=c(0.6, 1))
par(mar=c(0.8,0,0,0))
plot(crop(Grid100_2018,Grid1_2018[225,]),border="darkorange")
plot(Grid1_2018[225,],border="blue",add=T)
leg <- c(expression(paste("1 km"^"2")),expression(paste("0.01 km"^"2")))
legend("right" ,leg , col=c("blue","darkorange"),
       pch=c(0,0,0),bty="n",cex=1.65)
dev.off()

jpeg("./t3c.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
layout(matrix(c(1, 1,2, 2), nrow=2, byrow=TRUE),heights = c(2.2,0.35))
par(mar=c(0,0.8,0,0))
plot(zd_2018)
plot(Grid1_2018,border="blue",add=T)
plot(Grid100_2018,border="darkorange",add=T)
plot(Grid100_2018b,border="darkorange",add=T)
plot(infzone_2018,col="red",add=T)
plot(Grid1_2018,border="blue",add=T)
axis(side=1,at=seq(727000,780000,by=3000), pos=4270700)
axis(side=2,at=seq(4271000,4309000,by=3000), pos= 726500)
addscalebar(pos = "bottomright", plotepsg= 32630, padin=c(0.6, 1))
par(mar=c(0.8,0,0,0))
plot(crop(Grid100_2018b,Grid1_2018[225,]),border="darkorange")
plot(Grid1_2018[225,],border="blue",add=T)
leg <- c(expression(paste("1 km"^"2")),expression(paste("0.01 km"^"2")))
legend("right" ,leg , col=c("blue","darkorange"),
       pch=c(0,0,0),bty="n",cex=1.65)
dev.off()

# Table 4

# 3 phase design optimum values

# Grid 1 km2 
load(file="./results/three_phase/data1_50.RData") # Restriction 50%
load(file="./results/three_phase/data1_25.RData") # Restriction 25%
load(file="./results/three_phase/data1_15.RData") # Restriction 15%
load(file="./results/three_phase/data1_5.RData")  # Restriction 5%

replicates <- 100

aux1_50 <- list()
aux1_25 <- list()
aux1_15 <- list()
aux1_5 <- list()

celpos1_50 <- list()
celpos1_25 <- list()
celpos1_15 <- list()
celpos1_5 <- list()


for(i in 1:replicates){
  
  aux1_50[[i]] <- as.vector(data1_50[[i]]$nsamples)
  aux1_25[[i]] <- as.vector(data1_25[[i]]$nsamples)
  aux1_15[[i]] <- as.vector(data1_15[[i]]$nsamples)
  aux1_5[[i]] <- as.vector(data1_5[[i]]$nsamples)
  
  celpos1_50[[i]] <- sum(ifelse(data1_50[[i]]$positives >=1, 1, 0))
  celpos1_25[[i]] <- sum(ifelse(data1_25[[i]]$positives >=1, 1, 0))
  celpos1_15[[i]] <- sum(ifelse(data1_15[[i]]$positives >=1, 1, 0))
  celpos1_5[[i]] <- sum(ifelse(data1_5[[i]]$positives >=1, 1, 0))
  
}

n1_50 <- max(unlist(aux1_50))
n1_25 <- max(unlist(aux1_25))
n1_15 <- max(unlist(aux1_15))
n1_5 <- max(unlist(aux1_5))

c1_50 <- median(unlist(celpos1_50))
c1_25 <- median(unlist(celpos1_25))
c1_15 <- median(unlist(celpos1_15))
c1_5 <- median(unlist(celpos1_5))





# Grid 0.25 km2 

load(file="./results/three_phase/data500_50.RData") # Restriction 50%
load(file="./results/three_phase/data500_25.RData") # Restriction 25%
load(file="./results/three_phase/data500_15.RData") # Restriction 15%
load(file="./results/three_phase/data500_5.RData")  # Restriction 5%

aux500_50 <- list()
aux500_25 <- list()
aux500_15 <- list()
aux500_5 <- list()

celpos500_50 <- list()
celpos500_25 <- list()
celpos500_15 <- list()
celpos500_5 <- list()

for(i in 1:replicates){
  
  aux500_50[[i]] <- as.vector(data500_50[[i]]$nsamples)
  aux500_25[[i]] <- as.vector(data500_25[[i]]$nsamples)
  aux500_15[[i]] <- as.vector(data500_15[[i]]$nsamples)
  aux500_5[[i]] <- as.vector(data500_5[[i]]$nsamples)
  
  celpos500_50[[i]] <- sum(ifelse(data500_50[[i]]$positives >=1, 1, 0))
  celpos500_25[[i]] <- sum(ifelse(data500_25[[i]]$positives >=1, 1, 0))
  celpos500_15[[i]] <- sum(ifelse(data500_15[[i]]$positives >=1, 1, 0))
  celpos500_5[[i]] <- sum(ifelse(data500_5[[i]]$positives >=1, 1, 0))
}

n500_50 <- max(unlist(aux500_50))
n500_25 <- max(unlist(aux500_25))
n500_15 <- max(unlist(aux500_15))
n500_5 <- max(unlist(aux500_5))

c500_50 <- median(unlist(celpos500_50))
c500_25 <- median(unlist(celpos500_25))
c500_15 <- median(unlist(celpos500_15))
c500_5 <- median(unlist(celpos500_5))


# Grid 0.01 km2

load(file="./results/three_phase/data100_50.RData") # Restriction 50%
load(file="./results/three_phase/data100_25.RData") # Restriction 25%
load(file="./results/three_phase/data100_15.RData") # Restriction 15%
load(file="./results/three_phase/data100_5.RData")  # Restriction 5%

aux100_50 <- list()
aux100_25 <- list()
aux100_15 <- list()
aux100_5 <- list()

celpos100_50 <- list()
celpos100_25 <- list()
celpos100_15 <- list()
celpos100_5 <- list()

for(i in 1:replicates){
  
  aux100_50[[i]] <- as.vector(data100_50[[i]]$nsamples)
  aux100_25[[i]] <- as.vector(data100_25[[i]]$nsamples)
  aux100_15[[i]] <- as.vector(data100_15[[i]]$nsamples)
  aux100_5[[i]] <- as.vector(data100_5[[i]]$nsamples)
  
  celpos100_50[[i]] <- sum(ifelse(data100_50[[i]]$positives >=1, 1, 0))
  celpos100_25[[i]] <- sum(ifelse(data100_25[[i]]$positives >=1, 1, 0))
  celpos100_15[[i]] <- sum(ifelse(data100_15[[i]]$positives >=1, 1, 0))
  celpos100_5[[i]] <- sum(ifelse(data100_5[[i]]$positives >=1, 1, 0))
  

  }
  


n100_50 <- max(unlist(aux100_50))
n100_25 <- max(unlist(aux100_25))
n100_15 <- max(unlist(aux100_15))
n100_5 <- max(unlist(aux100_5))

c100_50 <- median(unlist(celpos100_50))
c100_25 <- median(unlist(celpos100_25))
c100_15 <- median(unlist(celpos100_15))
c100_5 <- median(unlist(celpos100_5))

phase1 <-  matrix(c(c1_50, n1_50, c1_25,n1_25, c1_15, n1_15, c1_5, n1_5,"-","-"), nrow=5, ncol=2, byrow=TRUE)

phase2 <-  matrix(c(c500_50, n500_50, c500_25, n500_25, c500_15, n500_15, c500_5, n500_5,"-","-"), nrow=5, ncol=2, byrow=TRUE)

phase3 <-  matrix(c(c100_50, n100_50, c100_25, n100_25, c100_15, n100_15, c100_5, n100_5,"-","-"), nrow=5, ncol=2, byrow=TRUE)

table4a <- rbind(phase1,phase2,phase3)
rownames(table4a) <- c("","j=1","","","","","j=0.25","","","","", "j=0.01","","","")
colnames(table4a) <- c("C*j,+","nj")




# 2 phase design optimum values


# Grid 1 km2 
load(file="./results/two_phase/data1_50.RData") # Restriction 50%
load(file="./results/two_phase/data1_25.RData") # Restriction 25%
load(file="./results/two_phase/data1_15.RData") # Restriction 15%
load(file="./results/two_phase/data1_5.RData")  # Restriction 5%

replicates <- 100

aux1_50 <- list()
aux1_25 <- list()
aux1_15 <- list()
aux1_5 <- list()

celpos1_50 <- list()
celpos1_25 <- list()
celpos1_15 <- list()
celpos1_5 <- list()


for(i in 1:replicates){
  
  aux1_50[[i]] <- as.vector(data1_50[[i]]$nsamples)
  aux1_25[[i]] <- as.vector(data1_25[[i]]$nsamples)
  aux1_15[[i]] <- as.vector(data1_15[[i]]$nsamples)
  aux1_5[[i]] <- as.vector(data1_5[[i]]$nsamples)
  
  celpos1_50[[i]] <- sum(ifelse(data1_50[[i]]$positives >=1, 1, 0))
  celpos1_25[[i]] <- sum(ifelse(data1_25[[i]]$positives >=1, 1, 0))
  celpos1_15[[i]] <- sum(ifelse(data1_15[[i]]$positives >=1, 1, 0))
  celpos1_5[[i]] <- sum(ifelse(data1_5[[i]]$positives >=1, 1, 0))
  
}

n1_50 <- max(unlist(aux1_50))
n1_25 <- max(unlist(aux1_25))
n1_15 <- max(unlist(aux1_15))
n1_5 <- max(unlist(aux1_5))

c1_50 <- median(unlist(celpos1_50))
c1_25 <- median(unlist(celpos1_25))
c1_15 <- median(unlist(celpos1_15))
c1_5 <- median(unlist(celpos1_5))


# Grid 0.01 km2

load(file="./results/two_phase/data100_50.RData") # Restriction 50%
load(file="./results/two_phase/data100_25.RData") # Restriction 25%
load(file="./results/two_phase/data100_15.RData") # Restriction 15%
load(file="./results/two_phase/data100_5.RData")  # Restriction 5%

aux100_50 <- list()
aux100_25 <- list()
aux100_15 <- list()
aux100_5 <- list()

celpos100_50 <- list()
celpos100_25 <- list()
celpos100_15 <- list()
celpos100_5 <- list()

for(i in 1:replicates){
  
  aux100_50[[i]] <- as.vector(data100_50[[i]]$nsamples)
  aux100_25[[i]] <- as.vector(data100_25[[i]]$nsamples)
  aux100_15[[i]] <- as.vector(data100_15[[i]]$nsamples)
  aux100_5[[i]] <- as.vector(data100_5[[i]]$nsamples)
  
  celpos100_50[[i]] <- sum(ifelse(data100_50[[i]]$positives >=1, 1, 0))
  celpos100_25[[i]] <- sum(ifelse(data100_25[[i]]$positives >=1, 1, 0))
  celpos100_15[[i]] <- sum(ifelse(data100_15[[i]]$positives >=1, 1, 0))
  celpos100_5[[i]] <- sum(ifelse(data100_5[[i]]$positives >=1, 1, 0))
  
  
}



n100_50 <- max(unlist(aux100_50))
n100_25 <- max(unlist(aux100_25))
n100_15 <- max(unlist(aux100_15))
n100_5 <- max(unlist(aux100_5))

c100_50 <- median(unlist(celpos100_50))
c100_25 <- median(unlist(celpos100_25))
c100_15 <- median(unlist(celpos100_15))
c100_5 <- median(unlist(celpos100_5))

phase1 <-  matrix(c(c1_50, n1_50, c1_25,n1_25, c1_15, n1_15, c1_5, n1_5,"-","-"), nrow=5, ncol=2, byrow=TRUE)

phase2 <-  matrix(c(c100_50, n100_50, c100_25, n100_25, c100_15, n100_15, c100_5, n100_5), nrow=4, ncol=2, byrow=TRUE)


table4b <- rbind(phase1,phase2)
rownames(table4b) <- c("","j=1","","","","","j=0.01","","")
colnames(table4b) <- c("C*j,+","nj")

table4 <-rbind(table4a,table4b)
table4

# Table 5

# Three phase design

load(file="./results/three_phase/data1_50.RData")
load(file="./results/three_phase/data500_50.RData")
load(file="./results/three_phase/data100_50.RData")

load("./data/binomial/data1.RData")
load("./data/binomial/data500.RData")
load("./data/binomial/data100.RData")


aux1_50 <- list()
aux500_50 <- list()
aux100_50 <- list()

for(i in 1:replicates){
  
  aux1_50[[i]] <- as.vector(data1_50[[i]]$nsamples)
  aux500_50[[i]] <- as.vector(data500_50[[i]]$nsamples)
  aux100_50[[i]] <- as.vector(data100_50[[i]]$nsamples)
}

n1_50 <- max(unlist(aux1_50))
n500_50 <- max(unlist(aux500_50))
n100_50 <- max(unlist(aux100_50))


phase_3 <- matrix(c(nrow(data1), n1_50, nrow(data1)* n1_50,
                    sum(ifelse(data1$positives >= 1 , 1, 0))*4, n500_50, sum(ifelse(data1$positives >= 1 , 1, 0))*4* n500_50,
                    sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))*25, n100_50, sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))*25*n100_50 ), 
                    nrow=3, ncol=3, byrow=TRUE)
# Two phase design

load(file="./results/two_phase/data1_50.RData")
load(file="./results/two_phase/data100_50.RData")

load("./data/binomial/data1.RData")


aux1_50 <- list()
aux100_50 <- list()

for(i in 1:replicates){
  
  aux1_50[[i]] <- as.vector(data1_50[[i]]$nsamples)
  aux100_50[[i]] <- as.vector(data100_50[[i]]$nsamples)
}

n1_50 <- max(unlist(aux1_50))
n100_50 <- max(unlist(aux100_50))


phase_2 <- matrix(c(nrow(data1), n1_50, nrow(data1)* n1_50,
                    sum(ifelse(data1$positives >= 1, 1, 0))*100, n100_50, sum(ifelse(data1$positives >= 1 , 1, 0))*100* n100_50),
                   nrow=2, ncol=3, byrow=TRUE)

# Current strategy 

Grid1_2018 <- readOGR("./data/grids", "Grid1_2018") # 1km x 1km grid (whole survey area)

zone100_2018 <- readOGR("./data/boundaries", "zone100_2018") # 100m x 100m grid boundaries

current <- matrix(c(round(((gArea(Grid1_2018)/10000)-(gArea(zone100_2018)/10000))/100), n1_50, round(((gArea(Grid1_2018)/10000)-(gArea(zone100_2018)/10000))/100) * n1_50,
                  round(gArea(zone100_2018)/10000), n100_50, round((gArea(zone100_2018)/10000))* n100_50),
                  nrow=2, ncol=3, byrow=TRUE)



table5 <- rbind(phase_3, phase_2, current)
rownames(table5) <- c("j=1","j=0.25","j=0.01","j=1","j=0.01","j=1","j=0.01")
colnames(table5) <- c("Cj","nj","Nj")
table5

# 3.2 Modelling the distribution of Xf incidence

# Table 6

# Model selected (based on WAIC)  resp ~ 1 + V

load(file="./results/incidence/models/mod13.Rdata")



table6a <- rbind(c(inla.zmarginal(mod13$marginals.fixed$`(Intercept)`)[c(1,2,5,3,7)],1-inla.pmarginal(0,mod13$marginals.fixed$`(Intercept)`)),
      
      c(inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),mod13$marginals.hyperpar$`Precision for V`))[c(1,2,5,3,7)], "-"))

# Model selected + climatic covariates resp ~ 1 + bio1 + bio7+ bio13 + V

load(file="./results/incidence/models/mod4.Rdata")


table6b <- rbind(c(inla.zmarginal(mod4$marginals.fixed$`(Intercept)`)[c(1,2,5,3,7)],1-inla.pmarginal(0,mod4$marginals.fixed$`(Intercept)`)),
      c(inla.zmarginal(mod4$marginals.fixed$bio1)[c(1,2,5,3,7)],1-inla.pmarginal(0,mod4$marginals.fixed$bio1)),
      c(inla.zmarginal(mod4$marginals.fixed$bio7)[c(1,2,5,3,7)],1-inla.pmarginal(0,mod4$marginals.fixed$bio7)),
      c(inla.zmarginal(mod4$marginals.fixed$bio13)[c(1,2,5,3,7)],1-inla.pmarginal(0,mod4$marginals.fixed$bio13)),
      c(inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),mod4$marginals.hyperpar$`Precision for V`))[c(1,2,5,3,7)], "-"))

table6 <- rbind(table6a, table6b)
rownames(table6) <- c("B0","sigma_v","B0","bio1","bio7","bio13","sigma_v")
colnames(table6) <- c("mean","sd","Q_0.5","Q_0.025","Q_0.975","P(·)>0")
table6


# Figure 4

# Figure 4a

jpeg("./4a.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
Grid1_2018$SPmean <- mod13$summary.random$V$mean
paleta<-get_col_regions()
par(mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=2, cex.lab=2, cex.main=2)
spplot(Grid1_2018, "SPmean", 
       col.regions=colorRampPalette(rev(paleta))(16),
       colorkey=list(space="bottom",height = 1), 
       par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"])
       ,at=seq(-2.15,6,by=0.6))

dev.off()

# Figure 4b

jpeg("./4b.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
Grid1_2018$presence_mean <- mod13$summary.fitted.values$mean
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018,"presence_mean", 
       col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=list(space="bottom",height = 1),
       par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),
       at=seq(0,0.34,by=0.02))
dev.off()

# Figure 4c

jpeg("./4c.jpeg",
     units="in", 
     width=8, 
     height=8, 
     pointsize=12, 
     res=300)
Grid1_2018$presence_sd <- mod13$summary.fitted.values$sd
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018,"presence_sd", 
       col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=list(space="bottom",height = 1),
       par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),
       at=seq(0,0.24,by=0.02))
dev.off()

# Effect of sampling intensity on Xf incidence estimates

# Database subsets models

load("./results/incidence/data_subsets_models/mod1_9.Rdata") # DS9
load("./results/incidence/data_subsets_models/mod1_23.Rdata") # DS23
load("./results/incidence/data_subsets_models/mod1_37.Rdata") # DS37
load("./results/incidence/data_subsets_models/mod1_51.Rdata") # DS51

# Parameters

replicates = 100
meanb0_9 = vector()
sdb0_9 = vector()
sigmav_9 =list()
meansigmav_9 = vector()
sdsigmav_9 = vector()

meanb0_23 = vector()
sdb0_23 = vector()
sigmav_23 = list()
meansigmav_23 = vector()
sdsigmav_23 = vector()

meanb0_37 = vector()
sdb0_37 = vector()
sigmav_37 = list()
meansigmav_37 = vector()
sdsigmav_37 = vector()

meanb0_51 = vector()
sdb0_51 = vector()
sigmav_51 = list()
meansigmav_51 = vector()
sdsigmav_51=vector()

for(i in 1:replicates){
  
  meanb0_9[i] = mod1_9[[i]]$summary.fixed$mean
  sdb0_9[i] = mod1_9[[i]]$summary.fixed$sd
  sigmav_9[[i]] <- inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),
                                                mod1_9[[i]]$marginals.hyperpar$`Precision for V`))
  
  meansigmav_9[i] = sigmav_9[[i]]$mean
  sdsigmav_9[i] = sigmav_9[[i]]$sd
  
  meanb0_23[i] = mod1_23[[i]]$summary.fixed$mean
  sdb0_23[i] = mod1_23[[i]]$summary.fixed$sd
  sigmav_23[[i]] <- inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),
                                                 mod1_23[[i]]$marginals.hyperpar$`Precision for V`))
  
  meansigmav_23[i] = sigmav_23[[i]]$mean
  sdsigmav_23[i] = sigmav_23[[i]]$sd
  
  
  meanb0_37[i] = mod1_37[[i]]$summary.fixed$mean
  sdb0_37[i] = mod1_37[[i]]$summary.fixed$sd
  sigmav_37[[i]] <- inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),
                                                  mod1_37[[i]]$marginals.hyperpar$`Precision for V`))
  
  meansigmav_37[i] = sigmav_37[[i]]$mean
  sdsigmav_37[i] = sigmav_37[[i]]$sd
  
  meanb0_51[i] = mod1_51[[i]]$summary.fixed$mean
  sdb0_51[i] = mod1_51[[i]]$summary.fixed$sd
  sigmav_51[[i]] <- inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),
                                                  mod1_51[[i]]$marginals.hyperpar$`Precision for V`))
  
  meansigmav_51[i] = sigmav_51[[i]]$mean
  sdsigmav_51[i] = sigmav_51[[i]]$sd
  
  
  
  
}


###########################

# Whole database model

load("./results/incidence_modeling/data_models/mod13.Rdata") 

# Parameters
b0 <- mod13$summary.fixed$mean
sv <- inla.zmarginal(inla.tmarginal(function(x) 1/sqrt(x),
                                    mod13$marginals.hyperpar$`Precision for V`))$mean



# Table 7


# DS9
beta0.ds9 <- c(mean(meanb0_9)-b0, mean(sdb0_9), sd(meanb0_9))
sigmav.ds9 <- c(mean(meansigmav_9)-sv, mean(sdsigmav_9), sd(meansigmav_9))

# DS23
beta0.ds23 <- c(mean(meanb0_23)-b0, mean(sdb0_23), sd(meanb0_23))
sigmav.ds23 <- c(mean(meansigmav_23)-sv, mean(sdsigmav_23), sd(meansigmav_23))

# DS37
beta0.ds37 <- c(mean(meanb0_37)-b0, mean(sdb0_37), sd(meanb0_37))
sigmav.ds37 <- c(mean(meansigmav_37)-sv, mean(sdsigmav_37), sd(meansigmav_37))

# DS51
beta0.ds51 <- c(mean(meanb0_51)-b0, mean(sdb0_51), sd(meanb0_51))
sigmav.ds51 <- c(mean(meansigmav_51)-sv, mean(sdsigmav_51), sd(meansigmav_51))


table7 <- cbind(rbind(beta0.ds9, sigmav.ds9, beta0.ds23, sigmav.ds23, beta0.ds37, sigmav.ds37, beta0.ds51, sigmav.ds51))
rownames(table7) <- c("B0","sigma_v","B0","sigma_v","B0","sigma_v","B0","sigma_v")
colnames(table7) <- c("Bias","SE","SD")
table7

# Figure 5

# Figure 5a

Grid1_2018$presence_mean <- mod13$summary.fitted.values$mean

# DS9
replicates=100
presence_mean=list()
for(i in 1:replicates){
  presence_mean[[i]] <- mod1_9[[i]]$summary.fitted.values$mean
}

Grid1_2018$presence_mean1 <- Reduce(`+`,presence_mean )/length(presence_mean) 
Grid1_2018$presence_meandif1 <- -Grid1_2018$presence_mean + Grid1_2018$presence_mean1
paleta <- colorRampPalette(c("blue","yellow","red"))(7)

jpeg(filename="./5a1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_meandif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.1,0.10,by=0.01))
#box(col = 'black')
dev.off()


# DS23
replicates=100
presence_mean=list()
for(i in 1:replicates){
  presence_mean[[i]] <- mod1_23[[i]]$summary.fitted.values$mean
}

Grid1_2018$presence_mean1 <- Reduce(`+`,presence_mean )/length(presence_mean) 
Grid1_2018$presence_meandif1 <- -Grid1_2018$presence_mean+Grid1_2018$presence_mean1

jpeg(filename="./5a2.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_meandif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.1,0.10,by=0.01))
#box(col = 'black')
dev.off()



# DS37
replicates=100
presence_mean=list()
for(i in 1:replicates){ 
  presence_mean[[i]] <- mod1_37[[i]]$summary.fitted.values$mean
}

Grid1_2018$presence_mean1 <- Reduce(`+`,presence_mean )/length(presence_mean) 
Grid1_2018$presence_meandif1 <- -Grid1_2018$presence_mean+Grid1_2018$presence_mean1

jpeg(filename="./5a3.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_meandif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.1,0.10,by=0.01))
#box(col = 'black')
dev.off()

# DS51
replicates=100
presence_mean=list()
for(i in 1:replicates){
  presence_mean[[i]] <- mod1_51[[i]]$summary.fitted.values$mean
}

Grid1_2018$presence_mean1 <- Reduce(`+`,presence_mean )/length(presence_mean) 
Grid1_2018$presence_meandif1 <- -Grid1_2018$presence_mean + Grid1_2018$presence_mean1


jpeg(filename="./5a4.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_meandif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=list(space="bottom",height = 1),par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.1,0.10,by=0.01))
#box(col = 'black')
dev.off()

# Figure 5b

Grid1_2018$presence_sd <- mod13$summary.fitted.values$sd

# DS9
replicates=100
presence_sd=list()
for(i in 1:replicates){
  presence_sd[[i]] <- mod1_9[[i]]$summary.fitted.values$sd
}

Grid1_2018$presence_sd1 <- Reduce(`+`,presence_sd )/length(presence_sd) 
Grid1_2018$presence_sddif1 <- -Grid1_2018$presence_sd+Grid1_2018$presence_sd1
paleta <- colorRampPalette(c("blue","yellow","red"))(7)

jpeg(filename="./5b1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_sddif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.023,0.035,by=0.003))
#box(col = 'black')
dev.off()


# DS23
replicates=100
presence_sd=list()
for(i in 1:replicates){
  presence_sd[[i]] <- mod1_23[[i]]$summary.fitted.values$sd
}

Grid1_2018$presence_sd1 <- Reduce(`+`,presence_sd )/length(presence_sd) 
Grid1_2018$presence_sddif1<- -Grid1_2018$presence_sd+Grid1_2018$presence_sd1
paleta<-colorRampPalette(c("blue","yellow","red"))(7)

jpeg(filename="./5b2.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_sddif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.023,0.035,by=0.003))
#box(col = 'black')
dev.off()


# DS37
replicates=100
presence_sd=list()
for(i in 1:replicates){
  presence_sd[[i]] <- mod1_37[[i]]$summary.fitted.values$sd
}

Grid1_2018$presence_sd1 <- Reduce(`+`,presence_sd )/length(presence_sd) 
Grid1_2018$presence_sddif1 <- -Grid1_2018$presence_sd+Grid1_2018$presence_sd1
paleta<-colorRampPalette(c("blue","yellow","red"))(7)
jpeg(filename="./5b3.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_sddif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.023,0.035,by=0.003))
#box(col = 'black')
dev.off()

# DS51
replicates = 100
presence_sd = list()
for(i in 1:replicates){
  presence_sd[[i]] <- mod1_51[[i]]$summary.fitted.values$sd
}

Grid1_2018$presence_sd1 <- Reduce(`+`,presence_sd )/length(presence_sd) 
Grid1_2018$presence_sddif1 <- -Grid1_2018$presence_sd+Grid1_2018$presence_sd1
paleta<-colorRampPalette(c("blue","yellow","red"))(7)

jpeg(filename="./5b4.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("presence_sddif1"),col.regions=colorRampPalette(rev(paleta))(20),
       colorkey=list(space="bottom",height = 1),par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-0.023,0.035,by=0.003))
#box(col = 'black')
dev.off()


# Figure 6c

Grid1_2018$Smean <- mod13$summary.random$V$mean 

# DS9
replicates=100
Smean=list()
for(i in 1:replicates){
  Smean[[i]] <- mod1_9[[i]]$summary.random$V$mean
}

Grid1_2018$Smean1 <- Reduce(`+`,Smean )/length(Smean) 
Grid1_2018$Smeandif1 <- -Grid1_2018$Smean + Grid1_2018$Smean1

paleta<-colorRampPalette(c("blue","yellow","red"))(7)
#paleta<-colorRampPalette(c("red","darkblue"))(8)

jpeg(filename="./5c1.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Smeandif1"),col.regions=colorRampPalette(rev(paleta))(12),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-2.7,1.03,by=0.3))
#box(col = 'black')
dev.off()



# DS23
replicates=100
Smean=list()
for(i in 1:replicates){
  Smean[[i]] <- mod1_23[[i]]$summary.random$V$mean
}

Grid1_2018$Smean1 <- Reduce(`+`,Smean )/length(Smean) 
Grid1_2018$Smeandif1 <- -Grid1_2018$Smean + Grid1_2018$Smean1


jpeg(filename="./6c2.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Smeandif1"),col.regions=colorRampPalette(rev(paleta))(12),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-2.7,1.03,by=0.3))
#box(col = 'black')
dev.off()


# DS37
replicates=100
Smean=list()
for(i in 1:replicates){
  Smean[[i]] <- mod1_37[[i]]$summary.random$V$mean
}

Grid1_2018$Smean1 <- Reduce(`+`,Smean )/length(Smean) 
Grid1_2018$Smeandif1 <- -Grid1_2018$Smean + Grid1_2018$Smean1

jpeg(filename="./5c3.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Smeandif1"),col.regions=colorRampPalette(rev(paleta))(12),
       colorkey=F,
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-2.7,1.03,by=0.3))
#box(col = 'black')
dev.off()

# DS51
replicates=100
Smean=list()
for(i in 1:replicates){
  Smean[[i]] <- mod1_51[[i]]$summary.random$V$mean
}

Grid1_2018$Smean1 <- Reduce(`+`,Smean )/length(Smean) 
Grid1_2018$Smeandif1 <- -Grid1_2018$Smean + Grid1_2018$Smean1

jpeg(filename="./5c4.jpeg", 
    units="in", 
    width=8, 
    height=6.5, 
    pointsize=12, 
    res=300)
par(mar=c(0,0,0,0), oma=c(0,0,0,0),cex.axis=0.7, cex.lab=0.7, cex.main=0.8)
spplot(Grid1_2018, c("Smeandif1"),col.regions=colorRampPalette(rev(paleta))(12),
       colorkey=list(space="bottom",height = 1),par.settings=list(fontsize=list(text=20)),
       xlim = c(Grid1_2018@bbox["x", "min"],
                Grid1_2018@bbox["x", "max"]),
       ylim = c(Grid1_2018@bbox["y", "min"],
                Grid1_2018@bbox["y", "max"]),at=seq(-2.7,1.03,by=0.3))
#box(col = 'black')
dev.off()







