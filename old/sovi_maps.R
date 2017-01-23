"
Functions for reading in and preparing data, plotting SoVI simulation & drop one maps and histograms.

"

library(maptools)
library(RColorBrewer)
library(classInt)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(stringr)
library(rgeos)
library(rgdal)

######################
## Data Preparation ##
######################

#read in SoVI variable names
SoVI_vars<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)$Attribute

#read in State queries
State_queries<-read.csv(paste0(getwd(),"/sovi_uncertainty/States/State_queries.csv"),stringsAsFactors=FALSE)

getSubset<-function(path,reg.type="FEMA",reg.name){
  
  "
  Gets queries for a subregion (FEMA region or State).

  path (character): current working directory

  reg.type (Character) : Region of interest type: 'FEMA' (default) or 'State'

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)
  "  

  #load the list of FEMA regions, separated by state
  if (reg.type=="FEMA"){
    region<-read.csv(paste0(path,"FEMA_queries.csv"),stringsAsFactors=FALSE)
  }else{
    region<-read.csv(paste0(path,"State_queries.csv"),stringsAsFactors=FALSE)  
  }
  
  ##Get the query for the region specified
  
  #If a state is specfied, need to change name to row index (integer)
  if (reg.type=="State"){
    reg.name<-which(region$State==reg.name)
  }
  
  query<-region$Query[reg.name]
  
  #Split query into values 
  query<-strsplit(query,"|g",fixed=TRUE)[[1]]

  #Remove leading "g" from first state FIPS 
  query[1]<-substr(query[1],2,3)
  
  return(query)
}

getRegionByQuery<-function(path,reg.name,reg.type="FEMA",cont_US=FALSE){
  
  "

  Returns a spatial polygons dataframe of counties for a specified region (USA, FEMA, or State)
  with baseline SoVI and SoVI rank.  

  path (character): current working directory

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  cont_US : set TRUE if mapping USA. Used to preserve rank values for 
            sim histograms. 
  "
  
  SoVI.sim<-readShapePoly(paste0(path,"/Build Data/USA_Counties_500k.shp"))

  sovi.prefix<-"SoVI_" #For reading in SoVI values
  
  #For distinguishing FEMA region/state inputs
  #region, for loading queries & naming outputs
  if (reg.type=="FEMA"){ #FEMA region specified 
    
    subpath<-"/sovi_uncertainty/FEMA_regions/"
    
    region<-paste0(reg.type,reg.name)
    sovi.prefix<-paste0(sovi.prefix,"FEMA_") #For reading in SoVI values
    
  }else if (reg.type=="USA"){ #Full USA specified
    
    subpath<-"/sovi_uncertainty/USA/"

    
  }else if (reg.type=="State"){ #State specified
    
    subpath<-"/sovi_uncertainty/States/"  
    
  }

  ##Subsetting by state(s) - only run "getSubset" for "FEMA" & "State"

  
  if (reg.type!="USA"){
    
    #Get FIPS code(s) for subsetting by state(s)
    subset.state<-getSubset(path=paste0(path,subpath),reg.type,reg.name)
    
    #Subset region
    SoVI.sim<-subset(SoVI.sim,SoVI.sim@data$STATE %in% subset.state)
    
  }
  
  #Read in computed SoVI values for region
  SoVI_reg<-read.csv(paste0(path,subpath,sovi.prefix,reg.name,"_all.csv"),stringsAsFactors=FALSE)
  
  SoVI.sim<-SoVI.sim[order(SoVI.sim$geoFIPS),]

  #Overwrite the existing SoVI values 
  SoVI.sim@data$sovi<-as.numeric(SoVI_reg$sovi)
  

  #Rank SoVI scores by county 
  SoVI.sim$sovi.rank<-rank(-SoVI.sim@data$sovi)
  
  #If "USA" and cont_US specified, subset continental USA for mapping
  if (reg.type=="USA" & cont_US){
    
    #Convert "STATE" to character
    SoVI.sim$STATE<-as.character(SoVI.sim$STATE)
    
    #Subset
    SoVI.sim<-SoVI.sim[SoVI.sim$STATE!="02" & SoVI.sim$STATE!="15",]
    
  }
  
  return(SoVI.sim)
}

###############
## SoVI Maps ##
###############

mapSoVI<-function(inData,outName,compare.state=NULL,condense.breaks=FALSE,map.quantile=FALSE,return.map_df=FALSE,borders=TRUE){
  
  "
  Plots SoVI for a given region of interest, either by raw values or by quantiles. 
  
  inData (Spatial Polygons Dataframe): Input region of interest
  
  outName (Character): output variable name, for file input/output and titling
  
  compare.state (Character): if not NULL, clips counties to a reference state - should only specify for USA & FEMA regions
  
  map.quantile (Logical): if TRUE, maps SoVI ranks rather than SoVI values (this mimics the HVRI SoVI maps)
  
  condense.breaks (Logical): if map.quantile=TRUE, if TRUE, displays only top 20% and bottom 20% of SoVI values
  
  return.map_df (Logical): if map.quantile=TRUE, returns the map dataframe object rather than plotting 
  
  borders (Logical): if TRUE (default), plots map borders 
  "
  
  #Prepare the data for mapping 
  
  if (map.quantile){ #Map SoVI quantiles - top 20%/bottom 20%
    
    #Generate breaks
    diff.breaks<-classIntervals(-inData$sovi,5,style="quantile")[[2]]
    
    #Generate labels 
    labels<-c("Top 20%","60% - 80%","40% - 60%","20% - 40%","Bottom 20%") 
    
    #Condense breaks if specified
    if (condense.breaks){
      
      diff.breaks<-classIntervals(-inData$sovi,3,style="quantile")[[2]]
      
      labels<-c("Top","Medium","Bottom")
      
    }
    
    #Assign breaks
    inData$value.bin<-cut(-inData$sovi,breaks=diff.breaks,labels=labels,include.lowest=TRUE)   
    
  }
  
  inData.name=outName
  
  inData.melt<-melt(inData@data,id.vars="geoFIPS")
  
  #fortify the sim_contig spatial dataframe
  map_geom<-fortify(inData,region="geoFIPS")
  
  #Subset the SoVI values
  
  if (map.quantile){ #Get binned sovi values
    
    sovi.all<-subset(inData.melt,variable %in% levels(inData.melt$variable)[which(levels(inData.melt$variable)=="value.bin")])
    
  }else{ #Subset raw sovi values 
    
    sovi.all<-subset(inData.melt,variable %in% levels(inData.melt$variable)[which(levels(inData.melt$variable)=="sovi")])
    
    #convert value field to numeric
    sovi.all$value<-as.numeric(sovi.all$value)
  }
  
  #Merge inData map 
  inData_map_df <- merge(map_geom, sovi.all, by.x="id", by.y="geoFIPS")
  
  #Generate title 
  map.title<-paste("SoVI by County,",inData.name)
  
  if (!is.null(compare.state)){ #Clip larger area of interest to specified state 
    
    #get state FIPS
    state<-State_queries[State_queries$State==compare.state,]$Query
    
    #Remove 1st character "g" from state query 
    state<-substr(state,2,nchar(state))
    
    #Clip inData
    inData_map_df<-inData_map_df[substr(inData_map_df$id,2,3)==state,]
    
    #Clip "sim_geom" (the boundaries layer)
    map_geom<-map_geom[substr(map_geom$id,2,3)==state,]
  }
  
  
  
  if (map.quantile){
    
    if (return.map_df){ #Get the map dataframe
      
      return(inData_map_df)
      
    }else{ #Plot the map dataframe
      
      #Make sure break labels are in order
      inData_map_df$value<-factor(inData_map_df$value,levels=labels)
      
      #plot original SoVI map 
      sovi.map<-ggplot(aes(long,lat,group=group),data=inData_map_df) + 
        geom_polygon(data=inData_map_df, aes(fill=value)) +
        coord_equal()+
        scale_fill_brewer(type="div",palette="RdBu")+
        guides(fill=guide_legend(title="SoVI"))+
        ggtitle(map.title)
      
      if (borders==TRUE){ #Plot map borders if specified
        
        sovi.map<-sovi.map+geom_path(data=map_geom, colour = "gray40", size = .5,alpha=.2)
        
      }
      
      sovi.map
      
    }
    
  }else{ #Map raw SoVI values 
    
    #Get minimum & maximum values 
    min.lim<-min(inData_map_df$value)
    max.lim<-max(inData_map_df$value)  
    
    if (return.map_df){ #Get the map dataframe
      
      return(inData_map_df)
      
    }else{ #Plot the map dataframe
      
      #For title, clean up inData.name for FEMA regions
      inData.name<-str_replace(inData.name,"_"," ")
      
      #plot original SoVI map 
      sovi.map<-ggplot(aes(long,lat,group=group),data=inData_map_df) + 
        geom_polygon(data=inData_map_df, aes(fill=value)) +
        coord_equal()+
        scale_fill_gradient2(limits=c(min.lim,max.lim),low="blue",mid="white",high="red",name="SoVI")+
        ggtitle(map.title)
      
      if (borders==TRUE){ #Plot map borders if specified
        
        sovi.map<-sovi.map+geom_path(data=map_geom, colour = "gray40", size = .5,alpha=.2)
        
      }
      
      sovi.map
      
    }
    
  }
  
}

######################################
## Margin-of-Error Simulation Tests ##
######################################

simStats<-function(inData,path,reg.type="FEMA",reg.name){
  
  "
  Calculates a series of descriptive statistics for 100 margin of error simulations of each variable 
  for a given area of interest and writes the results to an input spatial polygons dataframe. 

     -Mean of Changes in Rank: Mean of counties' net changes in rank across 100 MOE simulations of variable
     -Standard Deviation of Changes in Rank: Standard Deviation of counties' net changes in rank across 100 MOE simulations of variable
     -Mean of Absolute Value Changes in Rank: Mean of counties' changes in rank (positive or negative) across 100 MOE simulations of variable
     -Sum of Absolute Value Changes in Rank: Sum of counties' changes in rank (positive or negative) across 100 MOE simulations of variable
     -Maximum of Absolute Value Changes in Rank: Maximum  change in rank by county (positive or negative) across 100 MOE simulations of variable
     -Maximum - Minimum Changes in Rank: Ranges of counties' changes in rank (positive or negative) across 100 MOE simulations of variable

  inData (Spatial Polygons Dataframe): Input region of interest

  path (character): current working directory

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  "
  
  #Change region prefix based on reg.type
  if (reg.type=="FEMA"){ #FEMA region specified
    
    subpath<-"/sovi_uncertainty/FEMA_regions/"
    reg.prefix<-paste0("FEMA_",as.character(reg.name),"_")
    
  }else if (reg.type=="USA"){ #Full USA specified
    
    subpath<-"/sovi_uncertainty/USA/"
    reg.prefix<-"USA_" #For reading in SoVI values    
    
  }else if (reg.type=="State"){ #State specified
    
    subpath<-"/sovi_uncertainty/States/"
    reg.prefix<-paste0(reg.name,"_")
    
  }
  
  #Define set of output operations 
  operations<-c("rowMeans","mean.abs","sum.abs","maxMin","rowSds")
  
  #Define functions for absolute value of changes in rank
  mean.abs<-function(x) rowMeans(abs(x))
  sum.abs<-function(x) rowSums(abs(x))
  max.abs<-function(x) rowMaxs(abs(x))
  maxMin<-function(x) rowMaxs(x)-rowMins(x)
    
  #Iterate through SoVI inputs and generate summary stats 
  for (operation in operations){
    
    for (var.name in SoVI_vars){
             
      #read simulation output for current variable 
      #These need to be read in as text then converted to numeric because the default (factor) messes up the ranking 
      var.sim<-read.csv(paste0(path,subpath,reg.prefix,var.name,"_100.csv"),colClasses="character")
      
      #generate data frame for storing ranked simulated SoVI scores 
      ranks.sim<-data.frame(1:length(var.sim[,1]))
      
      #get changes in SoVI rank from baseline SoVI and assign to "ranks.sim"
      for (k in 2:length(var.sim)){
        
        #Rank the current SoVI simulation, converted to numeric
        sim.ranked<-rank(-as.numeric(var.sim[,k]))
        
        cur.sim<-k-2 #offset to match simulation number
        
        #Get the difference between the simulation rank and the original SoVI rank
        ranks.sim[[paste0("rank.chg_sim",cur.sim)]]<-inData@data$sovi.rank-sim.ranked
      }
      
      #remove junk first column
      ranks.sim<-ranks.sim[,2:length(ranks.sim)]
      
      #convert ranks.sim to matrix for matrixStats operations
      ranks.sim<-as.matrix(ranks.sim)
      
      #Assign columns to inData relative to the current operation 
      inData@data[[paste0("rankChg_",var.name,".",operation)]]<-do.call(operation,list(ranks.sim))
    }
  }
  
  return(inData)
}

plotSimHistogram<-function(path,inData,outName,opr="rowMeans"){
  
  "
  Plots a histogram of SoVI simulation results ('simStats' output) for a region interest by variable
  (in descending order of contribution to SoVI):

    -Mean of Changes in Rank: Mean of counties' net changes in rank across 100 MOE simulations of variable
    -Standard Deviation of Changes in Rank: Standard Deviation of counties' net changes in rank across 100 MOE simulations of variable
    -Mean of Absolute Value Changes in Rank: Mean of counties' changes in rank (positive or negative) across 100 MOE simulations of variable
    -Sum of Absolute Value Changes in Rank: Sum of counties' changes in rank (positive or negative) across 100 MOE simulations of variable
    -Maximum - Minimum Changes in Rank: Ranges of counties' changes in rank (positive or negative) across 100 MOE simulations of variable

  path (character): current working directory

  inData (Spatial Polygons Dataframe): Input region of interest

  outName (character): output variable name, for file input/output and titling

  opr: sets the operation in the simulation results to plot (default rowMeans)
  "
  
  #Get inData name as string for naming plot - need to extract from variable assignment
  inData.name=outName
  
  #Copy inData.name for file input 
  reg.prefix<-inData.name
  
  #Get sub-path for pulling SoVI net contrib
  if (grepl("FEMA",inData.name)){ #inData is a FEMA region 
    subpath<-"sovi_uncertainty/FEMA_regions/"
    
  }else if (grepl("USA",inData.name)){ #inData is full USA
    subpath<-"sovi_uncertainty/USA/"
    reg.prefix="USA"
    
  }else{ #inData is state 
    subpath<-"sovi_uncertainty/States/"
  }
  
  attr_contrib<-read.csv(paste0(path,subpath,reg.prefix,"_attribute+contrib_all.csv"),stringsAsFactors=FALSE)
  
  #Read in full SoVI vars - assign human readable attribute names
  var.names<-read.csv(paste0(path,"SoVI_vars.csv"),stringsAsFactors=FALSE)
  
  #Reassign attr_contrib attribute names as human readable names
  attr_contrib$attribute<-var.names$Desc
  
  #Get starting and ending positions of results to read
  startPos<-which(names(inData@data)==paste0("rankChg_MEDAGE_ACS.",opr))
  endPos<-which(names(inData@data)==paste0("rankChg_POPDENS.",opr))
  
  #Assign human-readable names to results
  names(inData@data)[startPos:endPos]<-attr_contrib$attribute  
  
  #Get index value for reading attribute names in for loop 
  attr.index<-1
  
  #Shorten sim result names & add net contribution for plotting
  for (i in startPos:endPos){
                
    #Pull the attribute contribution from CSV
    net_contrib<-as.character(round(attr_contrib[attr_contrib$attribute==names(inData@data)[i],]$rot_sum,2))

    #Combine variable & net contrib
    names(inData@data)[i]<-paste0(names(inData@data)[i],"\nnet contrib=",net_contrib)
    
    ##Also assign current attribute net contrib - we will reorder factor levels based on these and
    #need inData names and attributes to match 
    
    attr_contrib[attr.index,]$attribute<-paste0(attr_contrib[attr.index,]$attribute,"\nnet contrib=",net_contrib)
    
    #Increment attribute index
    attr.index<-attr.index+1
    
  }
  
  #Melt the inData dataframe
  sovi.sim_melt<-melt(inData@data[,c(which(names(inData@data)=="geoFIPS"),startPos:endPos)],id.vars="geoFIPS")
   
  #Change operation name for title 
  if (opr=="rowMeans"){
    opr="Mean"
  }else if (opr=="mean.abs"){
    opr="Mean Absolute Value"
  }else if (opr=="sum.abs"){
    opr="Sum of Absolute Value"
  }else if (opr=="maxMin"){
    opr="Maximum - Minimum"
  }else if (opr=="rowSds"){
    opr="Standard Deviation"
  }
  
  #Reorder variables based on full SoVI rank 
  attr_contrib$rank<-rank(-abs(as.numeric(attr_contrib$rot_sum)))
  
  attrs.ranked<-reorder(attr_contrib$attribute,attr_contrib$rank)
  
  #Reorder factor levels 
  sovi.sim_melt$variable<-factor(sovi.sim_melt$variable,levels=levels(attrs.ranked))
  
  #For title, clean up inData.name for FEMA regions
  inData.name<-str_replace(inData.name,"_"," ")
  
  #plot row means
  sovi.sim.plot<-ggplot(sovi.sim_melt,aes(x=as.numeric(value)))+
    geom_histogram(colour="#1c9099",fill="#f6eff7")+
    #xlim(-50,50)+
    facet_wrap(~variable,nrow=4,ncol=7,scales="free_x")+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    ggtitle(paste(opr,"SoVI Rank Change for",inData.name,"Counties,\n100 Margin of Error Simulations of SoVI Variables"))
  sovi.sim.plot                      
  
}

mapMaxMin<-function(path,inData,outName,by.var=NULL,borders=TRUE){
  
  "
  Maps the range of changes in rank for 100 margin of error simulations of SoVI variable(s) for a given region of interest.

  path (character): current working directory

  inData (Spatial Polygons Dataframe): Input region of interest ('simStats' output)

  outName (Character): output variable name, for file input/output and titling

  by.var (Character): If NULL (default) plots all MOE simulation results for all SoVI variables. 
    if not NULL, plots MOE simulation results a single variable.

   borders (Logical): if TRUE (default), plots map borders 
 
  "
  
  ##Prepare the data for mapping 
  
  #Get inData name as string for naming plot - need to extract from variable assignment
  inData.name=outName
  
  #Copy inData.name for file input 
  reg.prefix<-inData.name
  
  #Get sub-path for pulling SoVI net contrib
  if (grepl("FEMA",inData.name)){ #inData is a FEMA region 
    subpath<-"sovi_uncertainty/FEMA_regions/"
    
  }else if (grepl("USA",inData.name)){ #inData is full USA
    subpath<-"sovi_uncertainty/USA/"
    reg.prefix="USA"
    
    #Make sure AK and HI are removed, for plotting contiguous US
    inData<-inData[inData$STATE!="02" & inData$STATE!="15",]
    
  }else{ #inData is state 
    subpath<-"sovi_uncertainty/States/"
  }  
  
  #fortify the sim_contig spatial dataframe
  sim_geom<-fortify(inData,region="geoFIPS")
  
  #Melt the sim_contig values
  inData.melt<-melt(inData@data,id.vars="geoFIPS")
  
  operation<-"maxMin"
  opr.title<-"Maximum - Minimum Difference in"
  pal<-"RdPu" 
  pal.type<-"seq"
  
  #Get variables for mapping
  
  if (is.null(by.var)){ #Get list of all variables
    
    startPos<-which(levels(inData.melt$variable)==paste0("rankChg_",SoVI_vars[1],".",operation))
    endPos<-which(levels(inData.melt$variable)==paste0("rankChg_",SoVI_vars[length(SoVI_vars)],".",operation))
    
  }else{ #Plot only "by.var"
    
    startPos<-which(levels(inData.melt$variable)==paste0("rankChg_",by.var,".",operation))
    endPos<-which(levels(inData.melt$variable)==paste0("rankChg_",by.var,".",operation))
    
  }
    
  #subset by current operation
  map_subset<-subset(inData.melt,variable %in% levels(inData.melt$variable)[startPos:endPos])
  
  #convert value to numeric
  map_subset$value<-as.numeric(map_subset$value)
  
  for (j in startPos:endPos){
    
    #Get drop attribute name
    attr_mod<-strsplit(levels(map_subset$variable)[j],"rankChg_")[[1]][2]
    attr_mod<-strsplit(attr_mod,paste0(".",operation))[[1]][1]
    
#     attr_contrib<-read.csv(paste0(path,subpath,reg.prefix,"_attribute+contrib_all.csv"),stringsAsFactors=FALSE)
    
    #Pull the attribute contribution from CSV
#     attr_contrib<-attr_contrib[attr_contrib$attribute==attr_mod,]$rot_sum
    
    #Subset the SoVI values
    attr_sim<-subset(map_subset,variable %in% levels(map_subset$variable)[j])
    
    #convert value field to numeric
    attr_sim$value<-as.numeric(attr_sim$value)
    
    ##Set breaks - for consistency's sake aim for 5 intervals (6 break values),
    #but in some cases this will be lower and the number of intervals will need to be adjusted...
    
    #set default breaks number
    breaks.num<-5
    
    #get number of unique values 
    value.unique<-length(unique(attr_sim$value))
    
    if (value.unique<6){
      breaks.num<-value.unique-1
    }
    
    #Another safeguard - if we only have two unique values, then just assign the two unique values to "diff.breaks"
    if (breaks.num>1){
      
      diff.breaks<-classIntervals(unique(attr_sim$value),breaks.num,"quantile")[[2]]
      
    }else if (breaks.num==1){
      
      diff.breaks<-unique(attr_sim$value)
      
    }else{ #Only one unique value (breaks.num==0), skip the map 
      
      print(paste("Change in SoVI rank values are uniform for",outName,attr_mod,"(value=",unique(attr_sim$value),"). Did not plot map."))
      next
      
    }
    
    #Write intervals to labels
    labels<-vector()
    
    for (k in 1:(length(diff.breaks)-1)){
      labels<-c(labels,
                paste0('[',toString(diff.breaks[k]),' - ',toString(diff.breaks[k+1]),']'))
    }
    
    attr_sim$value.bin<-cut(attr_sim$value,breaks=diff.breaks,labels=labels,include.lowest=TRUE) 
    
    #Merge inData_geom with inData_map_rankChg attributes 
    inData_map_df <- merge(sim_geom, attr_sim, by.x="id", by.y="geoFIPS")

    #Get color ramp 
    fills<-rev(brewer.pal(length(levels(inData_map_df$value.bin)),pal))

    #Reverse labels for plotting legend
    inData_map_df$value.bin<-factor(inData_map_df$value.bin,levels=rev(labels))

    #Change attr_mod to human readable variable name 
    vars.list<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)
    attr_mod<-vars.list[vars.list$Attribute==attr_mod,]$Desc
    
    #For title, clean up inData.name for FEMA regions
    inData.name<-str_replace(inData.name,"_"," ")

    #Generate title 
#     map.title<-paste0(opr.title," Change in SoVI Rank by County, ",inData.name,"\n(100 Simulations of ",attr_mod,
#                       "; net contribution to SoVI=",as.character(round(attr_contrib,2)),")")

    map.title<-paste0(opr.title," Change in SoVI Rank by County, ",inData.name,"\n(100 Simulations of ",attr_mod)


    #plot - keep min/max interval limits consistent
    #this is accomplished by setting "drop=FALSE"" in the "scale_fill_brewer()" layer
    map.out<-ggplot(aes(long,lat,group=group),data=inData_map_df) + 
      geom_polygon(data=inData_map_df,aes(fill=value.bin))+
      coord_equal()+
      scale_fill_manual(values=fills,drop=FALSE,name="Range of Changes\nin SoVI Rank")+
      ggtitle(map.title)

    if (borders==TRUE){ #Plot map borders if specified
      
      map.out<-map.out+geom_path(data=sim_geom, colour = "gray40", size = .5,alpha=.2)
      
    }

    plot(map.out)
  }

}

mapSumAbs<-function(path,inData,outName,setBreaksEqual=FALSE,by.var=NULL,compare.state=NULL){
  
  "
  setBreaksEqual: If TRUE, displays equal interval break values; otherwise classifies by quantiles 
  
  by.var: if not NULL, maps a single variable 
  
  compare.state: if not NULL, clips counties to a reference state - should only specify for USA & FEMA regions
  "
  
  ##Prepare the data for mapping 
  
  #enable gpclib if it isn't already 
  gpclibPermit()
  
  #   #Get inData name as string for naming plot - need to extract from variable assignment
  #   inData.name<-toString(substitute(inData))
  #   inData.name<-strsplit(inData.name,", ")[[1]][2]
  inData.name=outName
  
  #Copy inData.name for file input 
  reg.prefix<-inData.name
  
  #Get sub-path for pulling SoVI net contrib
  if (grepl("FEMA",inData.name)){ #inData is a FEMA region 
    subpath<-"FEMA_regions/"
    
  }else if (grepl("USA",inData.name)){ #inData is full USA
    subpath<-"USA_all/"
    reg.prefix="USA_allUSA"
    
  }else{ #inData is state 
    subpath<-"States/"
  }  
  
  #   if (!is.null(compare.state)){ #Clip larger area of interest to specified state 
  #     
  #     #get state FIPS
  #     state<-State_queries[State_queries$State==compare.state,]$Query
  #     
  #     #Remove 1st character "g" from state query 
  #     state<-substr(state,2,nchar(state))
  #     
  #     #Clip inData
  #     inData<-inData[inData$STATE==state,]
  #     
  #   }
  
  #fortify the sim_contig spatial dataframe
  sim_geom<-fortify(inData,region="geoFIPS")
  
  #Melt the sim_contig values
  inData.melt<-melt(inData@data,id.vars="geoFIPS")
  
  operation<-"sum.abs"
  opr.title<-"Sum of Absolute Value"
  brewer.pal<-"YlOrRd"
  pal.type<-"seq"
  
  
  if (is.null(by.var)){ #Iterate through & map full set of variables 
    
    startPos<-which(levels(inData.melt$variable)==paste0("rankChg_",SoVI_vars[1],".",operation))
    endPos<-which(levels(inData.melt$variable)==paste0("rankChg_",SoVI_vars[length(SoVI_vars)],".",operation))
    
  }else{ #Only map specified variable
    
    startPos<-which(levels(inData.melt$variable)==paste0("rankChg_",by.var,".",operation))
    endPos=startPos
    
  }
  
  #subset by current operation
  map_subset<-subset(inData.melt,variable %in% levels(inData.melt$variable)[startPos:endPos])
  
  #convert value to numeric
  map_subset$value<-as.numeric(map_subset$value)
  
  for (j in startPos:endPos){
    
    #Get sim attribute name
    attr_mod<-strsplit(levels(map_subset$variable)[j],"rankChg_")[[1]][2]
    attr_mod<-strsplit(attr_mod,paste0(".",operation))[[1]][1]
    
    ##Get sim attribute contribution to SoVI
    #Read in CSV
    #     
    #     #Need to set a "reg.prefix" variable to distinguish USA from FEMA regions during file input
    #     if (reg.type=="USA"){
    #       
    #       reg.prefix="USA_allUSA"
    #       
    #     }else{ #keep reg.prefix as inData.name for FEMA and States
    #       
    #       reg.prefix=inData.name
    #       
    #     }
    
    #     attr_contrib<-read.csv(paste0(path,subpath,reg.prefix,"_attribute+contrib_all.csv"),stringsAsFactors=FALSE)
    
    
    #Pull the attribute contribution from CSV
    #     attr_contrib<-attr_contrib[attr_contrib$attribute==attr_mod,]$rot_sum
    
    #Subset the SoVI values
    attr_sim<-subset(map_subset,variable %in% levels(map_subset$variable)[j])
    
    #convert value field to numeric
    attr_sim$value<-as.numeric(attr_sim$value)
    
    ##Set breaks
    
    if (setBreaksEqual){ #Equal interval breaks
      
      ##for consistency's sake aim for 5 intervals (6 break values),
      #but in some cases this will be lower and the number of intervals will need to be adjusted...
      
      #set default breaks number
      breaks.num<-5
      
      #get number of unique values 
      value.unique<-length(unique(attr_sim$value))
      
      if (value.unique<6){
        breaks.num<-value.unique-1
      }
      
      #Another safeguard - if we only have one break value, then just assign the two unique values to "diff.breaks"
      if (breaks.num>1){
        diff.breaks<-classIntervals(attr_sim$value,breaks.num,"equal")[[2]]
        
      }else{
        
        diff.breaks<-unique(attr_sim$value)
        
      }
      
      diff.breaks<-classIntervals(attr_sim$value,breaks.num,"equal")[[2]]
      
      labels<-vector()
      
      for (k in 1:(length(diff.breaks)-1)){
        labels<-c(labels,
                  paste0('[',toString(diff.breaks[k]),' - ',toString(diff.breaks[k+1]),']'))
      }      
      
    }else{ #Quartile breaks
      
      diff.breaks<-classIntervals(attr_sim$value,n=4,style="quantile")[[2]]
      
      #Declare defined set of labels
      labels<-c("[0% - 25%]","[25% - 50%]","[50% - 75%]","[75%-100%]")
      
    }
    
    
    
    #Bin differences in rank by "diff.breaks"
    #If negative interval values exist, "cut" will assign them in the wrong order - need reverse 
    
    attr_sim$value.bin<-cut(attr_sim$value,breaks=diff.breaks,labels=labels,include.lowest=TRUE) 
    
    #Reassign zero values to the lowest interval - "cut" leaves these out...
    #   attr_sim[attr_sim$value==0,]$value.bin<-labels[1]
    
    #Merge inData_geom with inData_map_rankChg attributes 
    inData_map_df <- merge(sim_geom, attr_sim, by.x="id", by.y="geoFIPS")
    
    if (!is.null(compare.state)){ #Clip larger area of interest to specified state 
      
      #get state FIPS
      state<-State_queries[State_queries$State==compare.state,]$Query
      
      #Remove 1st character "g" from state query 
      state<-substr(state,2,nchar(state))
      
      #Clip inData
      inData_map_df<-inData_map_df[substr(inData_map_df$id,2,3)==state,]
      
      #Clip "sim_geom" (the boundaries layer)
      sim_geom<-sim_geom[substr(sim_geom$id,2,3)==state,]
      
    }  
    
    #Generate title 
    #     map.title<-paste0(opr.title," Change in SoVI Rank by County, ",inData.name,"\n(100 Simulations of ",attr_mod,
    #                       "; net contribution to SoVI=",as.character(round(attr_contrib,2)),")")
    
    map.title<-paste0(opr.title," Change in SoVI Rank by County, ",inData.name,"\n(100 Simulations of ",attr_mod)
    
    #plot - keep min/max interval limits consistent
    #this is accomplished by setting "drop=FALSE"" in the "scale_fill_brewer()" layer
    map.out<-ggplot(aes(long,lat,group=group),data=inData_map_df) + 
      geom_polygon(data=inData_map_df,aes(fill=value.bin))+
      coord_equal()+
      scale_fill_brewer(palette=brewer.pal,type=pal.type,drop=FALSE)+
      geom_path(data=sim_geom, colour = "gray40", size = .5,alpha=.2)+
      theme(legend.title=element_blank())+
      #     theme(legend.position="bottom") +
      ggtitle(map.title)
    plot(map.out)
  }
  
  return(map.out)
  
}


####################
## Drop-One Tests ##
####################

prep4DropOne<-function(path,reg.type="FEMA",reg.name,cont_US=FALSE){
  
  "
  Generates a spatial polygons dataframe for a given region of interest, and for each SoVI variable appends:
    -SoVI scores with variable drop
    -Ranks of SoVI scores with variable drop
    -Changes in rank from baseline SoVI with variable drop

  path (character): current working directory

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  cont_US : default FALSE; set to TRUE to map only continental USA 
  "
  
  #Generate a fresh dataframe without the sim stuff
  inData<-getRegionByQuery(path,reg.name,reg.type)
  
  #Change region prefix based on reg.type
  if (reg.type=="FEMA"){ #FEMA region specified
    
    subpath<-"/sovi_uncertainty/FEMA_regions/"
    reg.prefix<-paste0("SoVI_FEMA_",as.character(reg.name),"_")
    
  }else if (reg.type=="State"){ #State specified
    
    subpath<-"/sovi_uncertainty/States/"
    reg.prefix<-paste0("SoVI_",reg.name,"_")
    
  }else{ #USA specified
    
    subpath<-"/sovi_uncertainty/USA/"
    reg.prefix<-"SoVI_USA_" #For reading in SoVI values
    
  }

  for (var in SoVI_vars){

    inData[[paste0("drop_",var)]]<-as.numeric(as.vector(read.csv(paste0(path,subpath,reg.prefix,"drop_",var,".csv"))$sovi))
    
  }
    
  
  ###############################################
  ## Changes in SoVI Rank with Attribute Drops ##
  ###############################################
  
  #Rank the full SoVI values. Sort descending
  
  inData$sovi.rank<-rank(-inData@data$sovi)
  
  #Rank the drop-attribute SoVI values 
  
  for (i in which(names(inData)=="drop_MEDAGE_ACS"):which(names(inData)=="drop_POPDENS")){
    rank_field<-paste0("rank.",names(inData[i]))
    inData[[rank_field]]<-rank(-inData@data[,i])
  }
  
  #get the change in rank
  
  for (i in which(names(inData)=="rank.drop_MEDAGE_ACS"):which(names(inData)=="rank.drop_POPDENS")){
    fieldname<-paste0(names(inData[i]),"_chg")
    inData@data[[fieldname]]<-inData@data[,"sovi.rank"]-inData@data[,i]
  }
  
  #If "USA" and cont_US specified, subset continental USA for mapping
  if (reg.type=="USA" & cont_US){
    
    #Convert "STATE" to character
    inData$STATE<-as.character(inData$STATE)
    
    #Subset
    inData<-inData[inData$STATE!="02" & inData$STATE!="15",]
    
  }
  
  return(inData)
  
}

mapDropOne<-function(path,reg.type="FEMA",inData,outName,by.var=NULL,borders=TRUE){
  
  "
  path (character): current working directory

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  inData (Spatial Polygons Dataframe): Input region of interest ('simStats' output)

  outName (Character): output variable name, for file input/output and titling

  by.var (Character): If NULL (default) plots all MOE simulation results for all SoVI variables. 
    if not NULL, plots MOE simulation results a single variable.

   borders (Logical): if TRUE (default), plots map borders 
  "
  
  inData.name=outName

  #Get sub-path for pulling SoVI net contrib
  if (grepl("FEMA",inData.name)){ #inData is a FEMA region 
    subpath<-"sovi_uncertainty/FEMA_regions/"
    
  }else if (grepl("USA",inData.name)){ #inData is full USA
    subpath<-"sovi_uncertainty/USA/"
    
    #Make sure AK and HI are removed, for plotting contiguous US
    inData<-inData[inData$STATE!="02" & inData$STATE!="15",]
    
  }else{ #inData is state 
    subpath<-"sovi_uncertainty/States/"
  }
  
  ##Prepare the data for mapping 
  
  #fortify the input spatial dataframe
  map_geom<-fortify(inData,region="geoFIPS")
  
  #Melt the input dataframe 
  inData.melt<-melt(inData@data,id.vars="geoFIPS")
  
  if (is.null(by.var)){ #Get list of all variables
    
    startPos<-which(levels(inData.melt$variable)==paste0("rank.drop_",SoVI_vars[1],"_chg"))
    endPos<-which(levels(inData.melt$variable)==paste0("rank.drop_",SoVI_vars[length(SoVI_vars)],"_chg"))
    
  }else{ #Plot only "by.var"
    
    startPos<-which(levels(inData.melt$variable)==paste0("rank.drop_",by.var,"_chg"))
    endPos<-startPos
    
  }
  
  #subset the rank changes from the melted data frame 
  inData_map_rankChg<-subset(inData.melt,variable %in% levels(inData.melt$variable)[startPos:endPos])
  
  #convert value to character
  inData_map_rankChg$value<-as.numeric(inData_map_rankChg$value)
    
  ##Batch-plot dropped attributes
  
  for (i in startPos:endPos){
    
    #Get drop attribute name
    drop_attr<-strsplit(levels(inData_map_rankChg$variable)[i],"rank.drop_")[[1]][2]
    drop_attr<-strsplit(drop_attr,"_chg")[[1]][1]
    
    ##Get sim attribute contribution to SoVI
    #Read in CSV
    
    #Need to set a "reg.prefix" variable to distinguish USA from FEMA regions during file input
    if (reg.type=="USA"){
      
      reg.prefix="USA"
      
    }else{ #keep reg.prefix as inData.name for FEMA and States
      
      reg.prefix=inData.name
      
    }
    
    attr_contrib<-read.csv(paste0(path,subpath,reg.prefix,"_attribute+contrib_all.csv"),stringsAsFactors=FALSE)
    
    #Pull the attribute contribution from CSV
    attr_contrib<-as.character(round(attr_contrib[attr_contrib$attribute==drop_attr,]$rot_sum,2))
    
    #Subset the SoVI values
    sovi_drop<-subset(inData_map_rankChg,variable %in% levels(inData_map_rankChg$variable)[i])
    
    #convert value field to numeric
    sovi_drop$value<-as.numeric(sovi_drop$value)
    
    sovi_drop$geoFIPS<-as.character(sovi_drop$geoFIPS)
    
    #Assign break values based on rank changes 
    diff.breaks<-getDiffBreaks(unique(sovi_drop))
    labels<-getLabels(diff.breaks)
    
    #Bin differences in rank by "diff.breaks"
    sovi_drop$value.bin<-cut(sovi_drop$value,breaks=diff.breaks,labels=labels)
    
    #Merge inData_geom with inData_map_rankChg attributes 
    inData_map_df <- merge(map_geom, sovi_drop, by.x="id", by.y="geoFIPS")

    #Change drop_attr to human readable variable name 
    vars.list<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)
    drop_attr<-vars.list[vars.list$Attribute==drop_attr,]$Desc
        
    #For title, clean up inData.name for FEMA regions
    inData.name<-str_replace(inData.name,"_"," ")
    
    #For title, clean up inData.name for FEMA regions
    inData.name<-str_replace(inData.name,"_"," ")
    
    #Generate title 
    map.title<-paste0("Change in SoVI Rank by County, ",
                      inData.name,
                      "\n(drop ",
                      drop_attr,
                      "; weighted contribution to SoVI=",attr_contrib,")")
    
    fills<-brewer.pal(7,"Spectral")

    #Reverse labels for legend 
    inData_map_df$value.bin<-factor(inData_map_df$value.bin,levels=rev(labels))
    
    #plot - keep min/max interval limits consistent
    #this is accomplished by setting "drop=FALSE"" in the "scale_fill_brewer()" layer
    map.out<-ggplot(aes(long,lat,group=group),data=inData_map_df) + 
      geom_polygon(data=inData_map_df,aes(fill=value.bin))+
      coord_equal()+
      scale_fill_manual(values=fills,drop=FALSE,name="SoVI Rank Change")+
      ggtitle(map.title)

    if (borders==TRUE){ #Plot map borders if specified
      
      map.out<-map.out+geom_path(data=map_geom, colour = "gray40", size = .5,alpha=.2)
  
    }

    plot(map.out)
    
  }
}

getDiffBreaks<-function(inData){
  
  "
  Helper function #1 for getting legend values.
  
  Creates a set of breaks for SoVI rank changes with attribute drops 
  For regions/states, set the min & max by rank change values.
  
  inData (Spatial Polygons Dataframe): Input spatial polygons data frame for region of interest.
  
  "
  
  counties<-length(inData$geoFIPS)
  
  diff.breaks<-classIntervals(c(-counties:counties),7,"equal")[[2]]
  
  return(diff.breaks)
  
}

getLabels<-function(diff.breaks){
  
  "
  Helper function for getting labels from breaks (getDiffBreaks output)
  
  Create legend labels based on break values. 
  "
  
  labels<-vector()  
  
  for (k in 1:(length(diff.breaks)-1)){
    labels<-c(labels,
              paste0('[',toString(round(diff.breaks[k],2)),' - ',toString(round(diff.breaks[k+1],2)),']'))
  }
  
  return(labels)
  
}

#####################
## SoVI Comparison ##
#####################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  "
  Multiplot - from R Cookbook 
  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  
  Multiple plot function
  
  ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  - cols:   Number of columns in layout
  - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  
  If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  then plot 1 will go in the upper left, 2 will go in the upper right, and
  3 will go all the way across the bottom.
  "  
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mapCompareVar<-function(path,reg.name,reg.type,sovi_var,by.type="sim",map.quantile=FALSE,condense.breaks=FALSE,return.map_df=FALSE,clip.state=NULL,borders=TRUE){
  
  "
  Compares SoVI maps for a specified region (USA, FEMA, State) by drop-one or simulation. 

  path (character): current working directory

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  sovi_var (character): variable for comparison

  by.type (Character): if 'sim' (default) plots a comparison of simulated SoVI to original SoVI based on the 
           simulation of 'sovi_var' with the greatest observed movement in rank change. 

           if 'drop': plots a comparison of original SoVI to SoVI with sovi_var omitted. 

  map.quantile (Logical): if TRUE, maps SoVI ranks rather than SoVI values 

  condense.breaks (Logical): if map.quantile=TRUE, if TRUE, displays only top 20% and bottom 20% of SoVI values

  return.map_df (Logical): if map.quantile=TRUE, returns the map dataframe object rather than plotting 

  clip.state (Logical): if not NULL, clips larger regions of interest (USA, FEMA regions) to a specified state

  borders (Logical): if TRUE (default), plots map borders 
 
  "
  
  #Change region prefix based on reg.type
  if (reg.type=="FEMA"){ #FEMA region specified
    
    subpath<-"sovi_uncertainty/FEMA_regions/"
    reg.prefix<-paste0("FEMA_",as.character(reg.name),"_")
    
  }else if (reg.type=="State"){ #State specified
    
    subpath<-"sovi_uncertainty/States/"
    reg.prefix<-paste0("SoVI_",reg.name,"_")
    
  }else{ #USA specified
    
    subpath<-"sovi_uncertainty/USA/"
    reg.prefix<-"USA_" #For reading in SoVI values
    
  }
  
  #Get region
  reg<-getRegionByQuery(path,reg.name,reg.type,cont_US=FALSE)  
  
  if (by.type=="sim"){
    
    reg.max<-getSimDrop(path,reg,reg.name=reg.name,reg.type=reg.type,sovi_var=sovi_var,by.type=by.type)

    
    if (return.map_df==FALSE){ #Run a plot comparison
      
      #Get human-readable variable names 
      vars.list<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)
      sovi_var<-vars.list[vars.list$Attribute==sovi_var,]$Desc      
      
      #Clip USA to contiguous USA
      if (reg.type=="USA"){
        
        reg<-reg[reg$STATE!="02" & reg$STATE!="15",]
        reg.max<-reg.max[reg.max$STATE!="02" & reg.max$STATE!="15",]
      }
      
      #Expand FEMA region name
      if (reg.type=="FEMA"){
        reg.name<-paste("FEMA",toString(reg.name))
      }
      
      #Map SoVI for reg.max
      reg.max.sovi<-mapSoVI(reg.max,outName=paste(reg.name,"\nMaximum Simulation of",sovi_var),map.quantile=map.quantile,
                            condense.breaks=condense.breaks,borders=borders)
      
      #Map original SoVI
      reg.sovi<-mapSoVI(reg,outName=reg.name,map.quantile=map.quantile,condense.breaks=condense.breaks,borders=borders)

      #Compare to original reg SoVI map 
      multiplot(reg.sovi,reg.max.sovi,cols=2)
      
    }else{ #Get the map dfs and consolidate 

      mapCompareBins(reg=reg,reg.compare=reg.max,reg.name=reg.name,sovi_var=sovi_var,by.type="sim",condense.breaks=condense.breaks,
                     reg.type=reg.type,clip.state=clip.state)
      
    }      
     
  }else if (by.type=="drop"){
    
    reg.drop<-getSimDrop(path,reg,reg.name=reg.name,reg.type=reg.type,sovi_var=sovi_var,by.type=by.type)
    
    if (return.map_df==FALSE){ #Run a plot comparison
      
      #Get human-readable variable names 
      vars.list<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)
      sovi_var<-vars.list[vars.list$Attribute==sovi_var,]$Desc      
      
      #Clip USA to contiguous USA
      if (reg.type=="USA"){
        
        reg<-reg[reg$STATE!="02" & reg$STATE!="15",]
        reg.drop<-reg.drop[reg.drop$STATE!="02" & reg.drop$STATE!="15",]
      }
      
      #Expand FEMA region name
      if (reg.type=="FEMA"){
        reg.name<-paste("FEMA",toString(reg.name))
      }
    
      #Generate plot for comparison 
      multiplot(mapSoVI(reg,reg.name,map.quantile=map.quantile,condense.breaks=condense.breaks,borders=borders),
                mapSoVI(reg.drop,paste(reg.name,"- drop",sovi_var),map.quantile=map.quantile,condense.breaks=condense.breaks,borders=borders),
                cols=2)  
      
    }else{ #Get the map dfs and consolidate 
      
      mapCompareBins(reg=reg,reg.compare=reg.drop,reg.name=reg.name,sovi_var=sovi_var,by.type="drop",condense.breaks=condense.breaks,
                     reg.type=reg.type,clip.state=clip.state)
      
    }
    
  }
  
}

mapCompareBins<-function(reg,reg.compare,reg.type,reg.name,sovi_var,by.type,condense.breaks=FALSE,clip.state=NULL,borders=TRUE){
  
  "
  Helper function for mapCompareVar. Plots a quantile comparison (single map) of change between baseline SoVI 
  and Drop-One/Simulation test result for a given region of interest. 

  path (character): current working directory

  reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  sovi_var (character): variable for comparison

  by.type (Character): if 'sim' (default) plots a comparison of simulated SoVI to original SoVI based on the 
           simulation of 'sovi_var' with the greatest observed movement in rank change. 

           if 'drop': plots a comparison of original SoVI to SoVI with sovi_var omitted. 

  map.quantile (Logical): if TRUE, maps SoVI ranks rather than SoVI values 

  condense.breaks (Logical): if map.quantile=TRUE, if TRUE, displays only top 20% and bottom 20% of SoVI values

  return.map_df (Logical): if map.quantile=TRUE, returns the map dataframe object rather than plotting 

  "
  
  reg.map<-mapSoVI(reg,reg.name,map.quantile=TRUE,condense.breaks=condense.breaks,return.map_df=TRUE)
  reg_compare.map<-mapSoVI(reg.compare,reg.name,map.quantile=TRUE,condense.breaks=condense.breaks,return.map_df=TRUE)
  
  ##Assign comparison category to reg.map
  
  if (condense.breaks==TRUE){ #Condense breaks to 5 bins 
    
    reg.map$value<-ifelse(reg.map$value=="Top",3,
                              ifelse(reg.map$value=="Medium",2,1))
    
    reg_compare.map$value<-ifelse(reg_compare.map$value=="Top",3,
                                      ifelse(reg_compare.map$value=="Medium",2,1))
    
    
    #Get bin difference between original and result 
    reg.map$value.compare<-paste(reg_compare.map$value-reg.map$value)
    
    reg.map$value.compare<-factor(reg.map$value.compare,
                                  levels=c("2","1","0","-1","-2"))
    
    #Get colors for plotting 
    cols<-c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0")    
    
  }else{ #Keep full quantile range (7 bins)
    
    #Assign numerical values to bins and subtract 
    reg.map$value<-ifelse(reg.map$value=="Top 20%",5,
                              ifelse(reg.map$value=="60% - 80%",4,
                                     ifelse(reg.map$value=="40% - 60%",3,
                                            ifelse(reg.map$value=="20% - 40%",2,1))))
        
    reg_compare.map$value<-ifelse(reg_compare.map$value=="Top 20%",5,
                                  ifelse(reg_compare.map$value=="60% - 80%",4,
                                      ifelse(reg_compare.map$value=="40% - 60%",3,
                                             ifelse(reg_compare.map$value=="20% - 40%",2,1))))
    
    #Get bin difference between original and result 
    reg.map$value.compare<-paste(reg_compare.map$value-reg.map$value)
    
    reg.map$value.compare<-factor(reg.map$value.compare,
                                  levels=c("4","3","2","1","0","-1","-2","-3","-4"))
    
    #Get colors for plotting 
    cols<-c("#b2182b","#d6604d","#f4a582","#fddbc7",
            "#f7f7f7",
            "#d1e5f0","#92c5de","#4393c3","#2166ac")     
    
  }
  
  ##Set up for plotting contiguous USA 
  if (reg.type=="USA"){
    
    #Remove AK & HI 
    #Clip reg.map
    reg.map<-reg.map[!grepl("g02",reg.map$group) & !grepl("g15",reg.map$group),]
    
    #Clip boundary layer (reg)
    reg<-reg[reg$STATE!="02" & reg$STATE!="15",]
    
  }
  
  ##If FEMA or USA and clip.state is enabled, subset to state
  if (!is.null(clip.state) & (reg.type=="USA" | reg.type=="FEMA")){
    
    #Clip reg.map
    reg.map<-reg.map[grepl(clip.state,reg.map$group),]
    
    #Clip boundary layer (reg)
    reg<-reg[reg$STATE==substr(clip.state,2,3),]
    
  }  
  
  #Change sovi_var to human readable variable name 
  vars.list<-read.csv(paste0(getwd(),"/SoVI_vars.csv"),stringsAsFactors=FALSE)
  sovi_var<-vars.list[vars.list$Attribute==sovi_var,]$Desc
  
  #For title, clean up inData.name for FEMA regions
  if (reg.type=="FEMA"){
    reg.name<-paste("FEMA",as.character(reg.name))
  }        
  
  if (by.type=="sim"){
    by.type<-"Maximum Simulation of"
  }
  
  
  ##Plot reg.map based on value.compare
  sovi.map<-ggplot(aes(long,lat,group=group),data=reg.map) + 
    geom_polygon(data=reg.map, aes(fill=value.compare)) +
    coord_equal()+
    scale_fill_manual(values=cols,drop=FALSE,name="SoVI Change\n(Quantiles)")+
    ggtitle(paste(reg.name,"-",by.type,sovi_var))

  if (borders==TRUE){ #Plot map borders if specified
    
    sovi.map<-sovi.map+geom_path(data=reg, colour = "gray40", size = .5,alpha=.2)
    
  } 

  sovi.map
  
  return(sovi.map)
  
}

getSimDrop<-function(path,reg,reg.name,reg.type,sovi_var,by.type){
  
  "
  Helper function for mapCompareVar.

  For Simulation results, gets SoVI values for the margin-of-error simulation of a given SoVI variable with the most observed changes in rank. 

  For Drop-One results, gets SoVI values for the calculation of SoVI with a given variable omitted. 

  path (character): current working directory

  reg: input data

 reg.name (integer/character) : Name of region for file input/output. For FEMA regions, specify an integer
    (i.e. FEMA 4->4)

  reg.type (Character) : Region of interest type: 'USA', 'FEMA' (default), or 'State'

  sovi_var (character): variable for comparison

  by.type (Character): if 'sim' (default) plots a comparison of simulated SoVI to original SoVI based on the 
           simulation of 'sovi_var' with the greatest observed movement in rank change. 

           if 'drop': plots a comparison of original SoVI to SoVI with sovi_var omitted. 
  "
  
  #Change region prefix based on reg.type
  if (reg.type=="FEMA"){ #FEMA region specified
    
    subpath<-"sovi_uncertainty/FEMA_regions/"
    reg.prefix<-paste0("FEMA_",as.character(reg.name),"_")
    
  }else if (reg.type=="State"){ #State specified
    
    subpath<-"sovi_uncertainty/States/"
    reg.prefix<-paste0("SoVI_",reg.name,"_")
    
  }else{ #USA specified
    
    subpath<-"sovi_uncertainty/USA/"
    reg.prefix<-"USA_" #For reading in SoVI values
    
  }  
  
  if (by.type=="sim"){
    
    #Get simulation result
    if (reg.type=="USA" | reg.type=="FEMA"){
      
      reg.compare<-read.csv(paste0(path,subpath,reg.prefix,sovi_var,"_100.csv"),stringsAsFactors=FALSE)
      
    }else{ #State
      
      reg.compare<-read.csv(paste0(path,subpath,reg.name,"_",sovi_var,"_100.csv"),stringsAsFactors=FALSE)
      
    }  
    
    ##Sum simulations by total number of movements in rank (absolute value)
    
    #Create an empty vector for storing rank moves 
    rank.move<-vector()
    
    for (i in 2:ncol(reg.compare)){
      
      ##Get column rank change from original SoVI
      rankChg<-reg$sovi.rank-rank(-reg.compare[,i])
      
      #Append total sim movement to rank.move
      rank.move<-c(rank.move,sum(abs(rankChg)))
      
    }
    
    #Get index value of maximum move
    #Offset by one because sims start at index 2
    max.move<-reg.compare[,which.max(rank.move)+1]
    
    ##Map SoVI at max.move 
    
    #Reassign reg.compare as a copy of reg
    reg.compare<-reg
    
    #Assign max.move simulation as reg.compare SoVI values
    reg.compare@data$sovi<-max.move
    
    #Reassign SoVI rank
    reg.compare$sovi.rank<-rank(-reg.compare@data$sovi)
    
    #Rename if FEMA region
    if (reg.type=="FEMA"){
      reg.name<-paste("FEMA",as.character(reg.name))
    }
    
  }else{ # type 'drop' 
 
    ##Run drop-one on reg
    reg<-prep4DropOne(path,reg.type,reg.name)
    
    #Rename FEMA region
    if (reg.type=="FEMA"){
      reg.name<-paste("FEMA",as.character(reg.name))
    }
    
    #Create a copy of reg for sovi_var
    reg.compare<-reg
    
    #Assign sovi_var drop & rank to "reg.compare"
    reg.compare$sovi<-reg.compare[[paste0("drop_",sovi_var)]]
    reg.compare$sovi.rank<-reg.compare[[paste0("rank.drop_",sovi_var)]]
      
      
  }
  
  return(reg.compare)
 
}

