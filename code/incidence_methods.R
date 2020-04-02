###################################################################################################################
#
#                                                                                     
#   Filename    :	      incidence_methods.R  (UTF-8)  												  
#                                                                                     
#   Project     :       bioRxiv preprint "Tracking the outbreak. 
#                       An optimized delimiting survey strategy for Xylella fastidiosa
#                       in Alicante, Spain"                                                             
#   Authors     :       E. Lázaro et al.                                                              
#   Date        :       01.03.2020
#   Purpose     :       Implement methods described in Section 2.3 Modeling the distribution of Xf incidence
#                       as described in bioRxiv preprint
#																				  
#                                                 
#####################################################################################################################


# Packages
library(raster)
library(rgdal)
library(INLA)
library(spdep)
library(igraph)


setwd("C:/Users/MICOLOGIA_JOAQUIN/Dropbox/Dropbox/Xf/R/biorxiv_reproducibleresearch")

# Data

# Survey grids

Grid1_2018 <- readOGR("./data/grids", "Grid1_2018") # 1km x 1km grid (whole survey area)

# Survey data

load("./data/survey_data/xf2018.RData") # Data of 2018 delimiting Xf survey campaign in Alicante region

# Binomial data model

load("./data/binomial/data1.RData")

data1 <- data.frame(data1,V=c(1:length(data1[,1])), U=c(1:length(data1[,1]))) # V= spatial effect; U= independent random effect

# Neighborhood criterion

coords <- coordinates(Grid1_2018)

ID <- 1:dim(coords)[1] 

col.nb.0.all <- dnearneigh(coords, 0, 2500, row.names=ID)

nb2INLA("./data/neighborhood/grid1x1.graph", col.nb.0.all)


grid_grafo25 <- inla.read.graph("./data/neighborhood/grid1x1.graph")

# Model selection based on WAIC

source("./code/Bdiclcpomodel_stack.R") # selection model criteria

resp <- data1$positives

covariates <- c(paste("bio", c(1,7,13), sep=""), "f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)","f(U, model='iid')")

ntrials <- data1$nsamples


model_selection <- Bdiclcpomodel_stack(resp=resp, covariates=covariates,
                                Ntrials=ntrials,
                                database=data1, n=20,
                                family="binomial",
                                control.predictor=list(compute=FALSE), 
                                control.compute = list( 
                                  config=TRUE, dic=TRUE, cpo=TRUE, waic=TRUE), 
                                verbose=FALSE)


save(model_selection, file="./results/incidence/model_selection.RData")


# Models inference (the first 16 models based on WAIC)

# Model 1. resp ~ 1 + bio1 + bio7 + V

f1 <- positives ~ 1 + bio1 + bio7 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)

mod1 <- inla (f1, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod1_lcpo <- inla.cpo(mod1, verbose=F) # LCPO calculus. 

save(mod1,file="./results/incidence/models/mod1.Rdata")
save(mod1_lcpo,file="./results/incidence/models/mod1_lcpo.RData")


# Model 2.  resp ~ 1 + bio1 + bio7 + V + U

f2 <- positives ~ 1 + bio1 + bio7 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')

mod2 <- inla (f2, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod2_lcpo <- inla.cpo(mod2, verbose=F) # LCPO calculus. 

save(mod2,file="./results/incidence/models/mod2.Rdata")
save(mod2_lcpo,file="./results/incidence/models/mod2_lcpo.RData")

# Model 3.  resp ~ 1 + bio1 + bio7 + bio13 + V + U

f3 <- positives ~ 1 + bio1 + bio7 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')

mod3 <- inla (f3, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod3_lcpo <- inla.cpo(mod3, verbose=F) # LCPO calculus. 

save(mod3,file="./results/incidence/models/mod3.Rdata")
save(mod3_lcpo,file="./results/incidence/models/mod3_lcpo.RData")

# Model 4.    resp ~ 1 + bio1 + bio7 + bio13 + V

f4 <- positives ~ 1 + bio1 + bio7 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) 

mod4 <- inla (f4, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = FALSE)

mod4_lcpo <- inla.cpo(mod4, verbose=F) # LCPO calculus. 

save(mod4,file="./results/incidence/models/mod4.Rdata")
save(mod4_lcpo,file="./results/incidence/models/mod4_lcpo.RData")


# Model 5.    resp ~ 1 + bio7 + V + U

f5 <- positives ~ 1 + bio7  + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')

mod5 <- inla (f5, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose =F)

mod5_lcpo <- inla.cpo(mod5, verbose=F) # LCPO calculus. 

save(mod5,file="./results/incidence/models/mod5.Rdata")
save(mod5_lcpo,file="./results/incidence/models/mod5_lcpo.RData")


# Model 6.   resp ~ 1 + bio7 + V

f6 <- positives ~ 1 + bio7  + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)

mod6 <- inla (f6, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod6_lcpo <- inla.cpo(mod6, verbose= F) # LCPO calculus. 

save(mod6,file="./results/incidence/models/mod6.Rdata")
save(mod6_lcpo,file="./results/incidence/models/mod6_lcpo.RData")


# Model 7.    resp ~ 1 + bio7 + bio13 + V + U

f7 <- positives ~ 1 + bio7 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')

mod7 <- inla (f7, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod7_lcpo <- inla.cpo(mod7, verbose= F) # LCPO calculus. 

save(mod7,file="./results/incidence/models/mod7.Rdata")
save(mod7_lcpo,file="./results/incidence/models/mod7_lcpo.RData")


# Model 8.    resp ~ 1 + bio7 + bio13 + V

f8 <- positives ~ 1 + bio7 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) 

mod8 <- inla (f8, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)
 
mod8_lcpo <- inla.cpo(mod8, verbose= F) # LCPO calculus. 

save(mod8,file="./results/incidence/models/mod8.Rdata")
save(mod8_lcpo,file="./results/incidence/models/mod8_lcpo.RData")


# Model 9.    resp ~ 1 + bio1 + V + U

f9 <- positives ~ 1 + bio1 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid') 

mod9 <- inla (f9, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod9_lcpo <- inla.cpo(mod9, verbose= F) # LCPO calculus. 

save(mod9,file="./results/incidence/models/mod9.Rdata")
save(mod9_lcpo,file="./results/incidence/models/mod9_lcpo.RData")


# Model 10.       resp ~ 1 + bio1 + V

f10 <- positives ~ 1 + bio1 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)  

mod10 <- inla (f10, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
              control.compute = list(waic=TRUE, cpo= T), verbose = F)

mod10_lcpo <- inla.cpo(mod10, verbose= F) # LCPO calculus. 

save(mod10,file="./results/incidence/models/mod10.Rdata")
save(mod10_lcpo,file="./results/incidence/models/mod10_lcpo.RData")


# Model 11.     resp ~ 1 + bio1 + bio13 + V + U

f11 <- positives ~ 1 + bio1 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)  + f(U, model='iid')  

mod11 <- inla (f11, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod11_lcpo <- inla.cpo(mod11, verbose= F) # LCPO calculus. 

save(mod11,file="./results/incidence/models/mod11.Rdata")
save(mod11_lcpo,file="./results/incidence/models/mod11_lcpo.RData")


# Model 12.      resp ~ 1 + bio1 + bio13 + V

f12 <- positives ~ 1 + bio1 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)  

mod12 <- inla (f12, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod12_lcpo <- inla.cpo(mod12, verbose= F) # LCPO calculus. 

save(mod12,file="./results/incidence/models/mod12.Rdata")
save(mod12_lcpo,file="./results/incidence/models/mod12_lcpo.RData")


# Model 13.       resp ~ 1 + V 

f13 <- positives ~ 1  + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)  

mod13 <- inla (f13, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod13_lcpo <- inla.cpo(mod13, verbose= F) # LCPO calculus. 

save(mod13,file="./results/incidence/models/mod13.Rdata")
save(mod13_lcpo,file="./results/incidence/models/mod13_lcpo.RData")

# Model 14.       resp ~ 1 + bio13 + V + U

f14 <- positives ~ 1 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE)  + f(U, model='iid')  

mod14 <- inla (f14, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod14_lcpo <- inla.cpo(mod14, verbose= F) # LCPO calculus. 

save(mod14,file="./results/incidence/models/mod14.Rdata")
save(mod14_lcpo,file="./results/incidence/models/mod14_lcpo.RData")


# Model 15.      resp ~ 1 + bio13 + V

f15 <- positives ~ 1 + bio13 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) 

mod15 <- inla (f15, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod15_lcpo <- inla.cpo(mod15, verbose= F) # LCPO calculus. 

save(mod15,file="./results/incidence/models/mod15.Rdata")
save(mod15_lcpo,file="./results/incidence/models/mod15_lcpo.RData")


# Model 16.      resp ~ 1 + V + U

f16 <- positives ~ 1 + f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid') 

mod16 <- inla (f15, family = "binomial", data = data1, Ntrials = nsamples, control.predictor = list(compute=T, link = 1), 
               control.compute = list(waic=TRUE, cpo=TRUE), verbose = F)

mod16_lcpo <- inla.cpo(mod16, verbose= F) # LCPO calculus. 

save(mod16,file="./results/incidence/models/mod16.Rdata")
save(mod16_lcpo,file="./results/incidence/models/mod16_lcpo.RData")

####################################################################################


# Effect of sampling intensity 

# Creation data subsets: DS9, DS23, DS37, DS51.

Grid1_2018 <- spTransform(Grid1_2018, CRS("+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

celposneg <- over(xf2018[which(xf2018$reslab=="positive"|xf2018$reslab=="negative"),], as(Grid1_2018, "SpatialPolygons"))

celposneg <- as.vector(celposneg)

samples <- xf2018[-which(is.na(celposneg)),]

samples$cell <- celposneg[!is.na(celposneg)]

  replicates <- 100
  df <- list()
  data_s <- list()
  # sample1_9 <- list() # DS9
   # sample1_23 <- list() # DS23
   # sample1_37 <- list() # DS37
   sample1_51 <- list() # DS51
  positives <- list()
  negatives <- list()
  nsamples <- list()
  celpos <- list()
  cel <- list()
  samples_s <- list()
  samples_f <- list()
  # data1_9 <- list() # DS9
  # data1_23 <- list() # DS23
   # data1_37 <- list() # DS37
   data1_51 <- list() # DS51
  
  # j <- 9 
   # j <- 23 
   # j <- 37
   j <- 51
  
  
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
    # sample1_9[[i]] <- rbind(samples_f,data_s[[i]])
     # sample1_23[[i]] <- rbind(samples_f,data_s[[i]])
     # sample1_37[[i]] <- rbind(samples_f,data_s[[i]])
     sample1_51[[i]] <- rbind(samples_f,data_s[[i]])
    
    # positives[[i]] <- poly.counts(sample1_9[[i]][sample1_9[[i]]$reslab=="positive",], Grid1_2018)
    # negatives[[i]] <- poly.counts(sample1_9[[i]][sample1_9[[i]]$reslab=="negative",], Grid1_2018)
     # positives[[i]] <- poly.counts(sample1_23[[i]][sample1_23[[i]]$reslab=="positive",], Grid1_2018)
     # negatives[[i]] <- poly.counts(sample1_23[[i]][sample1_23[[i]]$reslab=="negative",], Grid1_2018)
     # positives[[i]] <- poly.counts(sample1_37[[i]][sample1_37[[i]]$reslab=="positive",], Grid1_2018)
     # negatives[[i]] <- poly.counts(sample1_37[[i]][sample1_37[[i]]$reslab=="negative",], Grid1_2018)
     positives[[i]] <- poly.counts(sample1_51[[i]][sample1_51[[i]]$reslab=="positive",], Grid1_2018)
     negatives[[i]] <- poly.counts(sample1_51[[i]][sample1_51[[i]]$reslab=="negative",], Grid1_2018)
    names(positives[[i]]) <- 1:nrow(Grid1_2018)
    names(negatives[[i]]) <- 1:nrow(Grid1_2018)
    nsamples[[i]] <- positives[[i]] + negatives[[i]]
    
    Pt <- coordinates(Grid1_2018)
    # data1_9[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data1_23[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     # data1_37[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
     data1_51[[i]] <- data.frame(x = Pt[,1], y = Pt[,2], positives = positives[[i]], nsamples = nsamples[[i]])
  }
  
# DS9
# save(data1_9, file="./data/binomial_subsets/data1_9.RData") # Binomial databases DS9

# DS23
# save(data1_23, file="./data/binomial_subsets/data1_23.RData") # Binomial databases DS23

# DS37
 # save(data1_37, file="./data/binomial_subsets/data1_37.RData") # Binomial databases DS37

# DS51
 save(data1_51, file="./data/binomial_subsets/data1_51.RData") # Binomial databases DS51



# Inference data subsets

# Databases

load("./data/binomial_subsets/data1_9.RData") # Binomial databases DS9
load("./data/binomial_subsets/data1_23.RData") # Binomial databases DS23
load("./data/binomial_subsets/data1_37.RData") # Binomial databases DS37
load("./data/binomial_subsets/data1_51.RData") # Binomial databases DS51

# Neighbourhood criterion

grid_grafo25 <- inla.read.graph("./data/neighborhood/grid1x1.graph")

replicates = 100
resp=list()
total=list()
# f1_9=list()
# f1_23=list()
 # f1_37=list()
 f1_51=list()
# mod1_9=list()
 # mod1_23=list()
 # mod1_37=list()
 mod1_51=list()

for (i in 1:replicates){
  
  # data1_9[[i]] <- data.frame(data1_9[[i]], V=c(1:length(data1_9[[i]][,1])), U=c(1:length(data1_9[[i]][,1])))
   # data1_23[[i]] <- data.frame(data1_23[[i]], V=c(1:length(data1_23[[i]][,1])), U=c(1:length(data1_23[[i]][,1])))
   # data1_37[[i]] <- data.frame(data1_37[[i]], V=c(1:length(data1_37[[i]][,1])), U=c(1:length(data1_37[[i]][,1])))
   data1_51[[i]] <- data.frame(data1_51[[i]], V=c(1:length(data1_51[[i]][,1])), U=c(1:length(data1_51[[i]][,1])))
   
  
  # f1_9 <- positives~ 1+ f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')
  # mod1_9[[i]]<-inla(f1_9, family="binomial", data=data1_9[[i]], Ntrials = nsamples, 
  #                  control.compute=list(dic=FALSE,cpo=TRUE, waic=TRUE), verbose = FALSE)
  # 
   # f1_23 <- positives~ 1+ f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')
   # mod1_23[[i]]<-inla(f1_23, family="binomial", data=data1_23[[i]], Ntrials = nsamples, 
   #                   control.compute=list(dic=FALSE,cpo=TRUE, waic=TRUE), verbose = FALSE)
  # 
   # f1_37 <- positives~ 1+ f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')
   # mod1_37[[i]]<-inla(f1_37, family="binomial", data=data1_37[[i]], Ntrials = nsamples, 
   #                   control.compute=list(dic=FALSE,cpo=TRUE, waic=TRUE), verbose = FALSE)
  # 
   f1_51 <- positives~ 1+ f(V, model='besag',graph=grid_grafo25, scale.model =TRUE) + f(U, model='iid')
   mod1_51[[i]]<-inla(f1_51, family="binomial", data=data1_51[[i]], Ntrials = nsamples, 
                     control.compute=list(dic=FALSE,cpo=TRUE, waic=TRUE), verbose = FALSE)
  
}

# save(mod1_9,file="./results/incidence/data_subsets_models/mod1_9.Rdata")
 # save(mod1_23,file="./results/incidence/data_subsets_models/mod1_23.Rdata")
 # save(mod1_37,file="./results/incidence/data_subsets_models/mod1_37.Rdata")
 save(mod1_51,file="./results/incidence/data_subsets_models/mod1_51.Rdata")

