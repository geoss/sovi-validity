

import os

os.chdir("/Users/Seth/Dropbox/SOVI_SIMULATOR/")

import sovi_simulator_simultaneous as sovi_simulator
import pandas as pd

path=os.path.join("/Users/Seth/Dropbox/SOVI_SIMULATOR/",'Build Data')

input_all=sovi_simulator.getInputNames()


#########
## USA ##
#########

outPath = os.path.join("/Users/Seth/Dropbox/SOVI_SIMULATOR/",'sovi_uncertainty','USA')

US_All = pd.read_csv(os.path.join(path,'sovi_inputs.csv'))
US_All.index = US_All.Geo_FIPS


#Compute full SoVi and simulations

sovi_simulator.computeSoVI(US_All,input_all,outPath,region='USA',attr='all',\
output_weights=True,simulation=True,sims=100)

#Compute drop-one 

for j in range(len(input_all)):
    
    #Reset input_names with each iteration
    input_names=sovi_simulator.getInputNames()
    
    attr_drop=input_names[j][0] #Get name of dropped attribute
    print "Dropping",attr_drop,"and computing SoVI..."
    input_names.pop(j) #Remove attribute from inputs
    
    #Compute SoVI 
    
    sovi_simulator.computeSoVI(US_All,input_names,outPath,'USA',attr=("drop_"+attr_drop),\
    output_weights=True,simulation=False,sims=None)

##################    
## FEMA Regions ##
##################
    
outPath = os.path.join(os.getcwd(),'sovi_uncertainty','FEMA_regions')

   ##get FEMA subregion queries    
FEMA_subs=pd.DataFrame.from_csv(os.path.join(os.getcwd(),outPath,"FEMA_queries.csv"))

#Just a test - keep adding to this with each iteration to 
#make sure we end up with the right number of counties (should equal contig US)
county_count=0    

for i in range(len(FEMA_subs)):

    #Get subregion name for output
    subr=FEMA_subs.ix[i].name    
    
    #Subset FEMA subregion
    data=US_All[US_All['Geo_FIPS'].str.contains(FEMA_subs['Query'][i])]
    
    print "Subset",subr+"...\n"
    
    #Compute full SoVi and simulations - up the number of simulations to 10000
    
    sovi_simulator.computeSoVI(data,input_all,outPath,region=subr,attr='all',\
    output_weights=True,simulation=True,sims=10000)
    
"""
    #Compute drop-one 

    for j in range(len(input_all)):
        
        #Reset input_names with each iteration
        input_names=sovi_simulator.getInputNames()
        
        attr_drop=input_names[j][0] #Get name of dropped attribute
        print "Dropping",attr_drop,"and computing SoVI..."
        input_names.pop(j) #Remove attribute from inputs
        
        #Compute SoVI 
        
        sovi_simulator.computeSoVI(data,input_names,outPath,subr,attr=("drop_"+attr_drop),\
        output_weights=True,simulation=False,sims=None)

    #Increment county count 
    county_count+=len(data)
    print "Got",str(county_count),"counties.\n"   
        
        
############        
## States ##
############ 

outPath = os.path.join(os.getcwd(),'sovi_uncertainty','States')      

#get state queries
State_subs=pd.DataFrame.from_csv(os.path.join(os.getcwd(),outPath,"State_queries.csv"))

for i in range(len(State_subs)):
   
    #Get subregion name for output    
    st=State_subs.ix[i].State    
    
    #Subset FEMA subregion
    data=US_All[US_All['Geo_FIPS'].str.contains(State_subs['Query'][i])]
    
    print "Subset",st+"...\n"
    
    #Compute full SoVi and simulations
   
    sovi_simulator.computeSoVI(data,input_all,outPath,region=st,attr='all',\
    output_weights=True,simulation=True,sims=100)

    #Compute drop-one 

    for j in range(len(input_all)):
        
        #Reset input_names with each iteration
        input_names=sovi_simulator.getInputNames()
        
        attr_drop=input_names[j][0] #Get name of dropped attribute
        print "Dropping",attr_drop,"and computing SoVI..."
        input_names.pop(j) #Remove attribute from inputs
        
        #Compute SoVI 
        
        sovi_simulator.computeSoVI(data,input_names,outPath,st,attr=("drop_"+attr_drop),\
        output_weights=True,simulation=False,sims=None)
"""