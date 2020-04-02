##########################################################################################################
#
#                                                                                     
#   Filename    :	      three_phase_strategy_methods.R (UTF-8)   												  
#                                                                                     
#   Project     :       bioRxiv preprint "Tracking the outbreak. 
#                       An optimized delimiting survey strategy for Xylella fastidiosa
#                       in Alicante, Spain"                                                             
#   Authors     :       E. Lázaro et al.                                                              
#   Date        :       01.03.2020
#   Purpose     :       Implement methods described in Section 3.1 Evaluation of delimiting strategies
#                       3 phase design (sequential adaptive strategy)
#																				  
#                                                 
###########################################################################################################


# Packages
library(raster)
library(rgdal)
library(rgeos)
library(spdep)
library(GISTools)
library(rlist)

setwd("C:/Users/MICOLOGIA_JOAQUIN/Dropbox/Dropbox/Xf/R/biorxiv_reproducibleresearch")

# Sequential adaptive strategy-->3 phase design

# Phase 1

# Data

# Survey data

load("./data/survey_data/xf2018.RData")

# Survey grids

Grid1_2018 <- readOGR("./data/grids", "Grid1_2018") # 1km x 1km grid (whole survey area)

# Binomial data

load("./data/binomial/data1.RData")

# Preparatoy phase

Grid1_2018 <- spTransform(Grid1_2018, CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

celposneg <- over(xf2018[which(xf2018$reslab=="positive"|xf2018$reslab=="negative"),], as(Grid1_2018, "SpatialPolygons"))

celposneg <- as.vector(celposneg)

samples <- xf2018[-which(is.na(celposneg)),]

samples$cell <- celposneg[!is.na(celposneg)]


# Phase 1

j <- 20 # number of samples per cell 
celpos = 0 # number of positive cells


 while(median(unlist(celpos)) < sum(ifelse(data1$positives >= 1 , 1, 0))){ # Restriction 50%
  # while(quantile(unlist(celpos),0.75) < sum(ifelse(data1$positives >= 1 , 1, 0))){  # Restriction 25%
 # while(quantile(unlist(celpos),0.85) < sum(ifelse(data1$positives >= 1 , 1, 0))){ # Restriction 15%
 # while(quantile(unlist(celpos),0.95) < sum(ifelse(data1$positives >= 1 , 1, 0))){ # Restriction 5%  
  
 rm(replicates, df, data_s, sample1_50, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) 
   # rm(replicates, df, data_s, sample1_25, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction 25%
 # rm(replicates, df, data_s, sample1_15, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction 15% 
 # rm(replicates, df, data_s, sample1_5, positives, negatives, nsamples, celpos, cel, samples_s, samples_f)  # Restriction  5% 
 
  
  replicates <- 100
  df <- list()
  data_s <- list()
   sample1_50 <- list()
    # sample1_25 <- list() # Restriction 25%
   # sample1_15 <- list() # Restriction 15%
   # sample1_5 <- list()  # Restriction 5%  
  positives <- list()
  negatives <- list()
  nsamples <- list()
  celpos <- list()
  cel <- list()
  samples_s <- list()
  samples_f <- list()
   data1_50 <- list()
    # data1_25 <- list()
   # data1_15 <- list()
   # data1_5 <- list()
 
  j <- j+1 
  
  # cells in which number of samples > j 
  cel <- as.vector(which(data1$nsamples >= j))
  
  samples$cel <- as.numeric(samples$cell %in% cel)
  
  # resampling samples
  samples_s <- samples[samples$cel == 1, ]
  
  # no resampling samples
  samples_f <- samples[samples$cel == 0,] 
  
  for(i in 1:replicates){
    
    # random sampling
    df[[i]] <- lapply(split(samples_s, samples_s@data$cell), function(x) x[sample(nrow(x), j), ])
    data_s[[i]] <- do.call("rbind", df[[i]]) 
    
    # new databases
     sample1_50[[i]] <- rbind(samples_f,data_s[[i]])
      # sample1_25[[i]] <- rbind(samples_f,data_s[[i]])
     # sample1_15[[i]] <- rbind(samples_f,data_s[[i]])
     # sample1_5[[i]] <- rbind(samples_f,data_s[[i]])
    
     positives[[i]] <- poly.counts(sample1_50[[i]][sample1_50[[i]]$reslab=="positive",], Grid1_2018)
     negatives[[i]] <- poly.counts(sample1_50[[i]][sample1_50[[i]]$reslab=="negative",], Grid1_2018)
      # positives[[i]] <- poly.counts(sample1_25[[i]][sample1_25[[i]]$reslab=="positive",], Grid1_2018)
      # negatives[[i]] <- poly.counts(sample1_25[[i]][sample1_25[[i]]$reslab=="negative",], Grid1_2018)
     # positives[[i]] <- poly.counts(sample1_15[[i]][sample1_15[[i]]$reslab=="positive",], Grid1_2018)
     # negatives[[i]] <- poly.counts(sample1_15[[i]][sample1_15[[i]]$reslab=="negative",], Grid1_2018)
     # positives[[i]] <- poly.counts(sample1_5[[i]][sample1_5[[i]]$reslab=="positive",], Grid1_2018)
     # negatives[[i]] <- poly.counts(sample1_5[[i]][sample1_5[[i]]$reslab=="negative",], Grid1_2018)
    names(positives[[i]]) <- 1:nrow(Grid1_2018)
    names(negatives[[i]]) <- 1:nrow(Grid1_2018)
    nsamples[[i]] <- positives[[i]] + negatives[[i]]
    
    Pt <- coordinates(Grid1_2018)
     data1_50[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
      # data1_25[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data1_15[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data1_5[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
  }
  
  
  for(i in 1:replicates){
     celpos[[i]] <-sum(ifelse(data1_50[[i]]$positives >= 1, 1, 0))
      # celpos[[i]] <-sum(ifelse(data1_25[[i]]$positives >= 1, 1, 0))
     # celpos[[i]] <-sum(ifelse(data1_15[[i]]$positives >= 1, 1, 0))
     # celpos[[i]] <-sum(ifelse(data1_5[[i]]$positives >= 1, 1, 0))
  }
  print(j)
}

# Restriction 50%
 save(data1_50, file="./results/three_phase/data1_50.RData") # Binomial databases

#  Restriction 25%
 # save(data1_25, file="./results/three_phase/data1_25.RData") # Binomial databases


#  Restriction 15%
# save(data1_15, file="./results/three_phase/data1_15.RData") # Binomial databases


#  Restriction 5%
# save(data1_5, file="./results/three_phase/data1_5.RData") # Binomial databases


# Optimum sampling intensity

aux1_50 <- list()
aux1_25 <- list()
aux1_15 <- list()
aux1_5 <- list()

for(i in 1:replicates){
  
aux1_50[[i]] <- as.vector(data1_50[[i]]$nsamples)
aux1_25[[i]] <- as.vector(data1_25[[i]]$nsamples)
aux1_15[[i]] <- as.vector(data1_15[[i]]$nsamples)
aux1_5[[i]] <- as.vector(data1_5[[i]]$nsamples)

}

opt1_50 <- max(unlist(aux1_50))
opt1_25 <- max(unlist(aux1_25))
opt1_15 <- max(unlist(aux1_15))
opt1_5 <- max(unlist(aux1_5))



# Phase 2

# Data

# Survey data

load("./data/survey_data/xf2018.RData")

# Survey grids

Grid500_2018 <- readOGR("./data/grids", "Grid500_2018") # 500m x 500m grid (whole survey area)
Grid500_2018 <- crop(Grid500_2018, Grid1_2018)
# Binomial data

load("./data/binomial/data500.RData") # Binomial data grid 500 x 500

load("./results/three_phase/data1_50.RData") # Binomial data phase 1 (100 replicates)

celpos<-list()

replicates <- 100

for(i in 1:replicates){
  
  celpos[[i]] <- sum(ifelse(data1_50[[i]]$positives >= 1, 1, 0))
}

da <- which(unlist(celpos) == sum(ifelse(data1$positives >= 1 , 1, 0)))

sample(da,1)


# Preparatoy phase

Grid1_2018_phase2 <- Grid1_2018[(which(data1_50[[sample(da,1)]]$positives >= 1)),]

Grid500_2018_phase2 <- crop(Grid500_2018,Grid1_2018_phase2)

trueCentroids = gCentroid(Grid500_2018_phase2, byid=TRUE)

coordinates <- trueCentroids@coords

Grid500_2018_phase2@data$XETRS89 <- coordinates[,1]

Grid500_2018_phase2@data$YETRS89 <- coordinates[,2]

Grid500_2018_phase2 <- Grid500_2018_phase2[,-c(1)]


Grid500_2018_phase2 <- spTransform(Grid500_2018_phase2, CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

celposneg <- over(xf2018[which(xf2018$reslab=="positive"|xf2018$reslab=="negative"),], as(Grid500_2018_phase2, "SpatialPolygons"))

celposneg <- as.vector(celposneg)

samples <- xf2018[-which(is.na(celposneg)),]

samples$cell <- celposneg[!is.na(celposneg)]

tab <- table(samples$cell)

idcell <- as.numeric(as.numeric(names(tab)))

nscell <- as.vector(tab)


# Phase 2

j <- 6 # number of samples per cell 
celpos = 0 # number of positive cells

  # while(median(unlist(celpos)) <  sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))){ # Restriction 50%
    # while(quantile(unlist(celpos),0.75) < sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))){ # Restriction 25%
   # while(quantile(unlist(celpos),0.85) < sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))){ # Restriction 15%
   while(quantile(unlist(celpos),0.95) < sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0))){ # Restriction 5%  
  
   # rm(replicates,df, data_s, sample500_50, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) 
   # rm(replicates,df, data_s, sample500_25, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction 25%
   # rm(replicates,df, data_s, sample500_15, positives, negatives, nsamples, celpos, cel, samples_s, samples_f)# Restriction 15% 
   rm(replicates,df, data_s, sample500_5, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction  5% 
  
  
  replicates <- 100
  df <- list()
  data_s <- list()
   # sample500_50 <- list()
   # sample500_25 <- list() # Restriction 25%
   # sample500_15 <- list() # Restriction 15%
   sample500_5 <- list()  # Restriction 5%  
  positives <- list()
  negatives <- list()
  nsamples <- list()
  celpos <- list()
  cel <- list()
  samples_s <- list()
  samples_f <- list()
    # data500_50 <- list()
   # data500_25 <- list()
   # data500_15 <- list()
   data500_5 <- list()
  
  j <- j+1 
  
  # cells in which number of samples > j 
  
  samples$cel <- as.numeric(samples$cell %in% idcell[ifelse(nscell > j, 1, 0) == 1])
  
  # resampling samples
  samples_s <- samples[samples$cel == 1, ]
  
  # no resampling samples
  samples_f <- samples[samples$cel == 0,] 
  
  for(i in 1:replicates){
    
    # random sampling
    df[[i]] <- lapply(split(samples_s, samples_s@data$cell), function(x) x[sample(nrow(x), j), ])
    data_s[[i]] <- do.call("rbind", df[[i]]) 
    
    # new databases
     # sample500_50[[i]] <- rbind(samples_f,data_s[[i]])
     # sample500_25[[i]] <- rbind(samples_f,data_s[[i]])
     # sample500_15[[i]] <- rbind(samples_f,data_s[[i]])
     sample500_5[[i]] <- rbind(samples_f,data_s[[i]])
    
     # positives[[i]] <- poly.counts(sample500_50[[i]][sample500_50[[i]]$reslab=="positive",], Grid500_2018_phase2[idcell,])
     # negatives[[i]] <- poly.counts(sample500_50[[i]][sample500_50[[i]]$reslab=="negative",], Grid500_2018_phase2[idcell,])
     # positives[[i]] <- poly.counts(sample500_25[[i]][sample500_25[[i]]$reslab=="positive",], Grid500_2018_phase2[idcell,])
     # negatives[[i]] <- poly.counts(sample500_25[[i]][sample500_25[[i]]$reslab=="negative",], Grid500_2018_phase2[idcell,])
     # positives[[i]] <- poly.counts(sample500_15[[i]][sample500_15[[i]]$reslab=="positive",], Grid500_2018_phase2[idcell,])
     # negatives[[i]] <- poly.counts(sample500_15[[i]][sample500_15[[i]]$reslab=="negative",], Grid500_2018_phase2[idcell,])
     positives[[i]] <- poly.counts(sample500_5[[i]][sample500_5[[i]]$reslab=="positive",], Grid500_2018_phase2[idcell,])
     negatives[[i]] <- poly.counts(sample500_5[[i]][sample500_5[[i]]$reslab=="negative",], Grid500_2018_phase2[idcell,])
    
    names(positives[[i]]) <- 1:nrow(Grid500_2018_phase2[idcell,])
    names(negatives[[i]]) <- 1:nrow(Grid500_2018_phase2[idcell,])
    nsamples[[i]] <- positives[[i]] + negatives[[i]]
    
    Pt <- Grid500_2018_phase2[idcell,]
     # data500_50[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data500_25[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data500_15[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     data500_5[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
  }
  
  
  for(i in 1:replicates){
    # celpos[[i]] <-sum(ifelse(data500_50[[i]]$positives >= 1, 1, 0))
       # celpos[[i]] <-sum(ifelse(data500_25[[i]]$positives >= 1, 1, 0))
     # celpos[[i]] <-sum(ifelse(data500_15[[i]]$positives >= 1, 1, 0))
     celpos[[i]] <-sum(ifelse(data500_5[[i]]$positives >= 1, 1, 0))
  }
  print(j)
}

# Restriction 50%
 # save(data500_50, file="./results/three_phase/data500_50.RData") # Binomial databases

#  Restriction 25%
# save(data500_25, file="./results/three_phase/data500_25.RData") # Binomial databases


#  Restriction 15%
 # save(data500_15, file="./results/three_phase/data500_15.RData") # Binomial databases


#  Restriction 5%
 save(data500_5, file="./results/three_phase/data500_5.RData") # Binomial databases


# Optimum sampling intensity

aux500_50 <- list()
aux500_25 <- list()
aux500_15 <- list()
aux500_5 <- list()

for(i in 1:replicates){
  
  aux500_50[[i]] <- as.vector(data500_50[[i]]$nsamples)
  aux500_25[[i]] <- as.vector(data500_25[[i]]$nsamples)
  aux500_15[[i]] <- as.vector(data500_15[[i]]$nsamples)
  aux500_5[[i]] <- as.vector(data500_5[[i]]$nsamples)
  
}

opt500_50 <- max(unlist(aux500_50))
opt500_25 <- max(unlist(aux500_25))
opt500_15 <- max(unlist(aux500_15))
opt500_5 <- max(unlist(aux500_5))


# Phase 3

# Data

# Survey data

load("./data/survey_data/xf2018.RData")

# Survey grids

Grid100_2018 <- readOGR("./data/grids", "Grid100_2018") # 100m x 100m grid (whole survey area)
Grid100_2018 <- crop(Grid100_2018, Grid1_2018)

# Binomial data

load("./data/binomial/data100.RData") # Binomial data grid 100 x 100

load("./results/three_phase/data500_50.RData") # Binomial data phase 2 (100 replicates)

celpos <- list()

replicates <- 100

for(i in 1:replicates){
  
  celpos[[i]] <- sum(ifelse(data500_50[[i]]$positives >= 1, 1, 0))
}

da <- which(unlist(celpos) ==  sum(ifelse((data500$positives >= 1) &  (!is.na(data500$positives)) * 1, 1, 0)))

replicate500 <- sample(da,1)





# Preparatoy phase


# Grid500_2018_phase3 <- Grid1_2018[(which(data500_50[[replicate500]]$positives >= 1)),]
# 
# Grid100_2018_phase3 <- crop(Grid100_2018,Grid500_2018_phase3)

X<-data500_50[[replicate500]][which(data500_50[[replicate500]]$positives>=1),]$XETRS89
Y<-data500_50[[replicate500]][which(data500_50[[replicate500]]$positives>=1),]$YETRS89
coor<-cbind(X, Y)
suma<-list()
for(i in 1:dim(Grid500_2018)[1]){
  suma[i]<-sum(coordinates(Grid500_2018)[i,1]==coor[,1] &coordinates(Grid500_2018)[i,2]==coor[,2])
}


Grid500_2018_phase3<-Grid500_2018[which(suma==1),]##89 celdas en el grid 500x500 que es el nÃºmero de celdas positivas
Grid100_2018_phase3<-crop(Grid100_2018, Grid500_2018_phase3)#2225 celdas en el grid 100x100



trueCentroids = gCentroid(Grid100_2018_phase3, byid=TRUE)

coordinates <- trueCentroids@coords

Grid100_2018_phase3@data$XETRS89 <- coordinates[,1]

Grid100_2018_phase3@data$YETRS89 <- coordinates[,2]

Grid100_2018_phase3 <- Grid100_2018_phase3[,-c(1)]


Grid100_2018_phase3 <- spTransform(Grid100_2018_phase3, CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))


celposneg <- over(xf2018[which(xf2018$reslab=="positive"|xf2018$reslab=="negative"),], as(Grid100_2018_phase3, "SpatialPolygons"))



celposneg <- as.vector(celposneg)

samples <- xf2018[which(!is.na(celposneg)),]

samples$cell <- celposneg[!is.na(celposneg)]

tab <- table(samples$cell)

idcell <- as.numeric(as.numeric(names(tab)))

nscell <- as.vector(tab)


# Phase 3

j <- 5 # number of samples per cell 
celpos = 0 # number of positive cells

 # while(median(unlist(celpos)) <  sum(ifelse((data100$positives >= 1) &  (!is.na(data100$positives)) * 1, 1, 0))){ # Restriction 50%
     # while(quantile(unlist(celpos),0.75) < sum(ifelse((data100$positives >= 1) &  (!is.na(data100$positives)) * 1, 1, 0))){ # Restriction 25%
   # while(quantile(unlist(celpos),0.85) < sum(ifelse((data100$positives >= 1) &  (!is.na(data100$positives)) * 1, 1, 0))){ # Restriction 15%
    while(quantile(unlist(celpos),0.95) < sum(ifelse((data100$positives >= 1) &  (!is.na(data100$positives)) * 1, 1, 0))){ # Restriction 5%  
  
   # rm(replicates,df, data_s, sample100_50, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) 
    # rm(replicates, df, data_s, sample100_25, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction 25%
   # rm(replicates, df, data_s, sample100_15, positives, negatives, nsamples, celpos, cel, samples_s, samples_f) # Restriction 15% 
   rm(replicates, df, data_s, sample100_5, positives, negatives, nsamples, celpos, cel, samples_s, samples_f)  # Restriction  5% 
  
  
  replicates <- 100
  df <- list()
  data_s <- list()
   # sample100_50 <- list()
    # sample100_25 <- list() # Restriction 25%
   # sample100_15 <- list() # Restriction 15%
   sample100_5 <- list()  # Restriction 5%  
  positives <- list()
  negatives <- list()
  nsamples <- list()
  celpos <- list()
  cel <- list()
  samples_s <- list()
  samples_f <- list()
   # data100_50 <- list()
    # data100_25 <- list()
   # data100_15 <- list()
   data100_5 <- list()
  
  j <- j+1 
  
  # cells in which number of samples > j 
  
  samples$cel <- as.numeric(samples$cell %in% idcell[ifelse(nscell > j, 1, 0) == 1])
  
  # resampling samples
  samples_s <- samples[samples$cel == 1, ]
  
  # no resampling samples
  samples_f <- samples[samples$cel == 0,] 
  
  for(i in 1:replicates){
    
    # random sampling
    df[[i]] <- lapply(split(samples_s, samples_s@data$cell), function(x) x[sample(nrow(x), j), ])
    data_s[[i]] <- do.call("rbind", df[[i]]) 
    
    # new databases
     # sample100_50[[i]] <- rbind(samples_f,data_s[[i]])
      # sample100_25[[i]] <- rbind(samples_f,data_s[[i]])
     # sample100_15[[i]] <- rbind(samples_f,data_s[[i]])
     sample100_5[[i]] <- rbind(samples_f,data_s[[i]])
    
     # positives[[i]] <- poly.counts(sample100_50[[i]][sample100_50[[i]]$reslab=="positive",], Grid100_2018_phase3[idcell,])
     # negatives[[i]] <- poly.counts(sample100_50[[i]][sample100_50[[i]]$reslab=="negative",], Grid100_2018_phase3[idcell,])
    # positives[[i]] <- poly.counts(sample100_25[[i]][sample100_25[[i]]$reslab=="positive",], Grid100_2018_phase3[idcell,])
    #   negatives[[i]] <- poly.counts(sample100_25[[i]][sample100_25[[i]]$reslab=="negative",], Grid100_2018_phase3[idcell,])
     # positives[[i]] <- poly.counts(sample100_15[[i]][sample100_15[[i]]$reslab=="positive",], Grid100_2018_phase3[idcell,])
     # negatives[[i]] <- poly.counts(sample100_15[[i]][sample100_15[[i]]$reslab=="negative",], Grid100_2018_phase3[idcell,])
     positives[[i]] <- poly.counts(sample100_5[[i]][sample100_5[[i]]$reslab=="positive",], Grid100_2018_phase3[idcell,])
     negatives[[i]] <- poly.counts(sample100_5[[i]][sample100_5[[i]]$reslab=="negative",], Grid100_2018_phase3[idcell,])
    names(positives[[i]]) <- 1:nrow(Grid100_2018_phase3[idcell,])
    names(negatives[[i]]) <- 1:nrow(Grid100_2018_phase3[idcell,])
    nsamples[[i]] <- positives[[i]] + negatives[[i]]
    
    Pt <- Grid100_2018_phase3[idcell,]
     # data100_50[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
      # data100_25[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data100_15[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     data100_5[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
  }
  
  
  for(i in 1:replicates){
     # celpos[[i]] <-sum(ifelse(data100_50[[i]]$positives >= 1, 1, 0))
     # celpos[[i]] <-sum(ifelse(data100_25[[i]]$positives >= 1, 1, 0))
     # celpos[[i]] <-sum(ifelse(data100_15[[i]]$positives >= 1, 1, 0))
     celpos[[i]] <-sum(ifelse(data100_5[[i]]$positives >= 1, 1, 0))
  }
  print(j)
  
}

# Restriction 50%
 # save(data100_50, file="./results/three_phase/data100_50.RData") # Binomial databases


#  Restriction 25%
# save(data100_25, file="./results/three_phase/data100_25.RData") # Binomial databases


#  Restriction 15%
 # save(data100_15, file="./results/three_phase/data100_15.RData") # Binomial databases

#  Restriction 5%
 save(data100_5, file="./results/three_phase/data100_5.RData") # Binomial databases


# Optimum sampling intensity

aux100_50 <- list()
aux100_25 <- list()
aux100_15 <- list()
aux100_5 <- list()

for(i in 1:replicates){
  
  aux100_50[[i]] <- as.vector(data100_50[[i]]$nsamples)
  aux100_25[[i]] <- as.vector(data100_25[[i]]$nsamples)
  aux100_15[[i]] <- as.vector(data100_15[[i]]$nsamples)
  aux100_5[[i]] <- as.vector(data100_5[[i]]$nsamples)
  
}

opt100_50 <- max(unlist(aux100_50))
opt100_25 <- max(unlist(aux100_25))
opt100_15 <- max(unlist(aux100_15))
opt100_5 <- max(unlist(aux100_5))
