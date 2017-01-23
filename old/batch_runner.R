"
Knits batch.rmd to generate all SoVI maps, simulation result histograms, simulation result maps, and drop one maps for 
contiguous USA, all FEMA regions, and all States. 

To execute in RStudio, select Code->Run Regions->Run All. 

Outputs are stored together under 'sovi_uncertainty' in the directory created for each subregion as:
  -an HTML file (all images) 
  -as graphics (PNG) in the 'figures' subdirectory (individual images) 
"

library(knitr)

path<-paste0(getwd(),"/")

createDirectories<-function(){
  
  "
  Create directories for outputs 
  "
  
  #FEMA Regions
  
  for (i in 1:10){
    dir.create(paste0(path,"sovi_uncertainty/FEMA_regions/FEMA_",as.character(i)))
  }
  
  #States
  states<-read.csv(paste0(path,"sovi_uncertainty/States/State_queries.csv"),stringsAsFactors=FALSE)$State
  
  for (state in states){
    dir.create(paste0(path,"sovi_uncertainty/States/",state))
  }
  
}

batchUSA<-function(){
  
  reg.type="USA"
  reg.name="USA"  
  
  #Set working directory
  setwd(dir=paste0(path,"sovi_uncertainty/USA"))
  
  outName<-"USA" #For naming output figs
  
  knit2html(input=paste0(path,"batch.Rmd"),
            output=paste0(getwd(),"/output.html"))
  
}

batchAllFEMAState<-function(){
  
  "Generates outputs for all FEMA Regions & States using batchOutput
  "
  
  #Get list of states
  states<-read.csv(paste0(path,"sovi_uncertainty/States/State_queries.csv"),stringsAsFactors=FALSE)$State
  
  #Iterate through states & FEMA regions and generate outputs
  for (i in 1:10){
    batchOutput(i,states[i])
  }
  
}

batchFEMAState<-function(FEMA.name,State.name){

  #################
  ## FEMA Region ##
  #################
  
  reg.type="FEMA"
  reg.name=FEMA.name
  
  #Set working directory
  setwd(dir=paste0(path,"sovi_uncertainty/FEMA_regions/FEMA_",as.character(FEMA.name)))
  
  outName=paste0(reg.type,"_",as.character(reg.name)) #For naming output figs
  
  knit2html(input=paste0(path,"batch.Rmd"),
            output=paste0(getwd(),"/output.html"))
  
  ############
  ## States ##
  ############
  
  reg.type="State"
  reg.name=State.name
  
  #Set working directory
  setwd(dir=paste0(path,"sovi_uncertainty/States/",State.name))
  
  outName=reg.name #For naming output figs
  
  knit2html(input=paste0(path,"batch.Rmd"),
            output=paste0(getwd(),"/output.html"))
  
  
}


#########
## Run ##
#########

createDirectories()
batchUSA()
batchAllFEMAState()