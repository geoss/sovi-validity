source(paste0(getwd(),"/sovi_maps.R"))

path=getwd()

##Load regions 

#USA 
USA<-getRegionByQuery(path,reg.name="USA",reg.type="USA")
plot(USA)

#Contiguous USA (don't know whether to keep this...)
cont_USA<-getRegionByQuery(path,reg.name="USA",reg.type="USA",cont_US=TRUE)
plot(cont_USA)

#FEMA region
FEMA_8<-getRegionByQuery(path,reg.name=8,reg.type="FEMA")
plot(FEMA_8)

#State
SD<-getRegionByQuery(path,reg.name="South Dakota",reg.type="State")
plot(SD)

##############
## Map SoVI ##
##############

##USA

#Raw values
mapSoVI(cont_USA,"USA")
mapSoVI(cont_USA,"USA",borders=FALSE)

#By quantiles
mapSoVI(cont_USA,"USA",map.quantile = TRUE)

#By quantiles (condensed)
mapSoVI(cont_USA,"USA",map.quantile = TRUE,condense.breaks=TRUE)
mapSoVI(cont_USA,"USA",map.quantile = TRUE,condense.breaks=TRUE,borders=FALSE)

##FEMA

#Raw values
mapSoVI(FEMA_8,"FEMA 8")
mapSoVI(FEMA_8,"FEMA 8",borders=FALSE)

#By quantiles (condensed)
mapSoVI(FEMA_8,"FEMA 8",map.quantile = TRUE,condense.breaks=TRUE)
mapSoVI(FEMA_8,"FEMA 8",map.quantile = TRUE,condense.breaks=TRUE,borders=FALSE)

##State

#Raw values
mapSoVI(SD,"South Dakota")
mapSoVI(SD,"South Dakota",borders=FALSE)

#By quantiles (condensed)
mapSoVI(SD,"South Dakota",map.quantile = TRUE,condense.breaks=TRUE)
mapSoVI(SD,"South Dakota",map.quantile = TRUE,condense.breaks=TRUE,borders=FALSE)


######################
## Drop-one figures ##
######################

##USA

#Prepare data
USA<-prep4DropOne(path,reg.type="USA",reg.name="USA")

#Map drop-one (just one variable)
mapDropOne(path,reg.type="USA",inData=USA,outName="USA",by.var="QAGEDEP_ACS")

#Map drop-one (all variables)
mapDropOne(path,reg.type="USA",inData=USA,outName="USA")

##FEMA Region
FEMA_4<-prep4DropOne(path,"FEMA",4)

#Map drop-one (just one variable)
mapDropOne(path,reg.type="FEMA",inData=FEMA_4,outName="FEMA_4",by.var="QSERV_ALT")

#Map drop-one (all variables)
mapDropOne(path,reg.type="FEMA",inData=FEMA_4,outName="FEMA_4")

##State
Georgia<-prep4DropOne(path,"State","Georgia")

#Map drop-one (just one variable)
mapDropOne(path,reg.type="State",inData=Georgia,outName="Georgia",by.var="QNATAM_ACS")

#Map drop-one (all variables)
mapDropOne(path,reg.type="State",inData=Georgia,outName="Georgia")


########################
## Simulation figures ##
########################

##USA
USA<-simStats(getRegionByQuery(path,"USA","USA"),path,"USA","USA")

#Get histograms for variable movement (default - mean rank change by county)
plotSimHistogram(path,inData=USA,outName="USA")

#Map max-min rank change (just one variable)
mapMaxMin(path,inData=USA,outName="USA",by.var="QRICH200K")

#Map max-min rank change (all variables)
mapMaxMin(path,inData=USA,outName="USA")


##FEMA region
FEMA_6<-simStats(getRegionByQuery(path,6,"FEMA"),path,"FEMA",6)

#Get histograms for variable movement (default - mean rank change by county)
plotSimHistogram(path,inData=FEMA_6,outName="FEMA_6")

#Map max-min rank change (just one variable)
mapMaxMin(path,inData=FEMA_6,outName="FEMA_6",by.var="QNRRES_ACS")

#Map max-min rank change (all variables)
mapMaxMin(path,inData=FEMA_6,outName="FEMA_6")


##State
Texas<-simStats(getRegionByQuery(path,"Texas","State"),path,"State","Texas")

#Get histograms for variable movement (default - mean rank change by county)
plotSimHistogram(path,inData=Texas,outName="Texas")

#Map max-min rank change (just one variable)
mapMaxMin(path,inData=Texas,outName="Texas",by.var="QNRRES_ACS")

#Map max-min rank change (all variables)
mapMaxMin(path,inData=Texas,outName="Texas")


##########################
## Quantile Change Maps ##
##########################

### USA ###

##Simulation

#Full quantile breaks
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QRICH200K",by.type="sim",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)

#Condensed breaks
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QRICH200K",by.type="sim",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QRICH200K",by.type="sim",
              return.map_df=FALSE)

##Drop-One

#Full quantile breaks
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)


#Condensed breaks 
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name="USA",reg.type="USA",sovi_var="QESL_ALT",by.type="drop",
              return.map_df=FALSE)

### FEMA Region ###

#Full quantile breaks
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QNATAM_ACS",by.type="sim",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)

#Condensed breaks
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QNATAM_ACS",by.type="sim",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QNATAM_ACS",by.type="sim",
              return.map_df=FALSE)

##Drop

#Full quantile breaks
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)


#Condensed breaks 
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name=3,reg.type="FEMA",sovi_var="QESL_ALT",by.type="drop",
              return.map_df=FALSE)

### State ###

##Sim 

#Full quantile breaks
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QAGEDEP_ACS",by.type="sim",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)

#Condensed breaks
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QAGEDEP_ACS",by.type="sim",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QAGEDEP_ACS",by.type="sim",
              return.map_df=FALSE)

##Drop

#Full quantile breaks
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=FALSE,return.map_df=TRUE)


#Condensed breaks 
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QESL_ALT",by.type="drop",
              map.quantile=TRUE,condense.breaks=TRUE,return.map_df=TRUE)

#Side-by-side SoVI comparison
mapCompareVar(path,reg.name="Virginia",reg.type="State",sovi_var="QESL_ALT",by.type="drop",
              return.map_df=FALSE)







