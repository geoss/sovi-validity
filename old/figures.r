########################################################
## Script to Produce Figures in Variable-wise SoVI Paper
## Seth Spielman
## Adapted from Joe Tuccillo's "heatmaps.rmd"
## June 2015
########################################################

##Libraries
library(ggplot2)
library(reshape2)


#Set working path
setwd("~/Dropbox/SOVI_SIMULATOR")

path<-paste0(getwd(),"/")

#Get list of SoVI variables
SoVI_vars<-as.vector(read.csv(paste0(path,"paper/SoVI_vars.csv"))$Attribute)

#Get a List of FEMA regions
FEMA_regions<-as.vector(read.csv(paste0(path,"data/FEMA_regions/FEMA_queries.csv"))$Subregion)

#Add "US_All" to represent all US counties
FEMA_regions<-c("US_All", FEMA_regions)

#List of States
States<-c(read.csv(paste0(path,"data/States/State_queries.csv"),stringsAsFactors=FALSE)$Abbrev)

#matrices for storing rank and net contrib by variable by geography
weights_mat<-matrix(nrow=length(SoVI_vars),ncol=length(FEMA_regions))
rownames(weights.mat)<-SoVI_vars
colnames(weights.mat)<-FEMA_regions

ranks_mat<-matrix(nrow=length(SoVI_vars),ncol=length(FEMA_regions))
rownames(ranks.mat)<-SoVI_vars
colnames(ranks.mat)<-FEMA_regions


