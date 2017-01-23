###########################
##Script to create variable
###########################

import os as os
import numpy as np
import pandas as pd
from scipy.stats.mstats import zscore as ZSCORE
from scipy.stats import rankdata
from scipy.stats import spearmanr 
from code.spss_pca import SPSS_PCA


local_path = '/Users/Seth/'
#local_path = '/Users/dfolch/'

os.chdir(local_path+'Dropbox/SoVI_var_wise_paper/code')
path=local_path+'/Dropbox/SoVI_var_wise_paper'
outPath =local_path+'/Dropbox/SoVI_var_wise_paper/data'

US_All = pd.read_csv(path+'/data/input/sovi_inputs.csv')
US_All.index = US_All.Geo_FIPS


# attribute name and expected influence on vulnerability
input_names = [['MEDAGE_ACS','pos','person'],
               ['BLACK_ACS','pos','person'],   
               ['QNATAM_ACS','pos','person'],  
               ['QASIAN_ACS','pos','person'],  
               ['QHISP_ACS','pos','person'],   
               ['QAGEDEP_ACS','pos','person'], 
               ['QPUNIT_ACS','pos','person'],  
               ['PRENTER_ACS','pos','hu'], 
               ['QNRRES_ACS','pos','person'],  
               ['QFEMALE_ACS','pos','person'], 
               ['QFHH_ACS','pos','hu'],    
               ['QUNOCCHU_ACS','pos','hu'],
               ['PERCAP_ALT','neg','person'],  
               ['QESL_ALT','pos','person'],    
               ['QCVLUN','pos','person'], 
               ['QPOVTY','pos','person'],
               ['QMOHO','pos','hu'],
               ['QED12LES_ALT','pos','person'],
               ['QFEMLBR','pos','person'], 
               ['QEXTRCT_ALT','pos','person'], 
               ['QSERV_ALT','pos','person'],   
               ['QSSBEN','pos','hu'],
               ['QNOAUTO_ALT','pos','hu'], 
               ['QFAM','neg','person'], 
               ['QRICH200K','neg','hu'],
               ['MDGRENT_ALT','neg','hu'], 
               ['MHSEVAL_ALT','neg','hu'], 
               ['POPDENS','pos','person']] 

#Get attribute names
attr_names=[j[0] for j in input_names]
#cols = [c for c in US_All.columns if c.find('_SE') == -1]

attr_names.append('Geo_FIPS')
#US_All = US_All.dropna(axis=0) #two counties misisng data in state 15 and 48
US_All = US_All[attr_names]
US_All['stateID'] = US_All.Geo_FIPS.str.slice(0,3,1)
attr_names.remove('Geo_FIPS')

# input data prep
# --swap signs of the attributes expected to have a "negative" affect on vulnerability
# --take z-score of data
for name, sign, sample in input_names:
    if sign == 'neg':
        US_All[name] = -US_All[name].values
    elif sign == 'pos':
        pass
    else:
        raise Exception, "problem"


########################################
##SOVI FOR FEMA REGIONS
#########################################
#Build FEMA subRegions Dict values= state ID's
FEMA_subs= {'FEMA_1':['g23','g50','g33','g25','g09','g44']}
FEMA_subs['FEMA_2'] = ['g36','g34']
FEMA_subs['FEMA_3'] = ['g42','g10','g11','g24','g51','g54']
FEMA_subs['FEMA_4'] = ['g21','g47','g37','g28','g01','g13','g45','g12']
FEMA_subs['FEMA_5'] = ['g27','g55','g26','g17','g18','g39']
FEMA_subs['FEMA_6'] = ['g35','g48','g40','g05','g22']
FEMA_subs['FEMA_7'] = ['g31','g19','g20','g29']
FEMA_subs['FEMA_8'] = ['g30','g38','g56','g46','g49','g08']
FEMA_subs['FEMA_9'] = ['g06','g32','g04']
FEMA_subs['FEMA_10'] = ['g53','g41','g16']

#Dict to hold variable loadings
varContrib = {}

#Multiindexed DataFrame to hold all FEMA SOVI Scores
geoLevels = US_All.Geo_FIPS
femaLevels = FEMA_subs.keys()
geoLabels = []
femaLabels = []
for f in femaLevels:
    femaRegionIndexes = US_All[US_All['stateID'].isin(FEMA_subs[f])].index.values
    geoLabels.extend([US_All.index.get_loc(i) for i in femaRegionIndexes])
    femaLabels.extend(np.repeat(femaLevels.index(f), len(femaRegionIndexes)))

US_femaSub_Multi_Index = pd.MultiIndex(levels=[femaLevels, geoLevels], 
                                    labels=[femaLabels, geoLabels], 
                                    names=['FEMA_Region', 'Geo_FIPS'])

FEMA_Region_Sovi_Score = pd.DataFrame(index=US_femaSub_Multi_Index, columns=['sovi', 'rank']) 

for i in FEMA_subs:
    
    #Subset FEMA subregion
    FEMARegionData=US_All[US_All['stateID'].isin(FEMA_subs[i])]

    # compute SoVI
    inputData = FEMARegionData.drop(['Geo_FIPS','stateID'], axis = 1, inplace = False)
    pca = SPSS_PCA(inputData, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=FEMARegionData.Geo_FIPS, columns=['sovi'])
    attrib_contribution = pca.weights_rot.sum(1)
    
    FEMA_Region_Sovi_Score.loc[i, 'sovi'] = sovi_actual.values
    #ADD RANKabs(sovi_actual).apply(rankdata, axis=0, method='average')
        
    ##Write attribute contribution output     
    #Generate dictionary for all net loadings by variable and region
    varContrib[i]=zip(attr_names,attrib_contribution.tolist())

#######################
##Compute National SoVI
#######################
# compute SoVI
inputData = US_All.drop(['Geo_FIPS','stateID'], axis = 1, inplace = False)
pca = SPSS_PCA(inputData, reduce=True, varimax=True)
sovi_actual = pca.scores_rot.sum(1)
sovi_actual = pd.DataFrame(sovi_actual, index=US_All.Geo_FIPS, columns=['sovi'])
attrib_contribution = pca.weights_rot.sum(1)
US_All_Full_Sovi_Rank = abs(sovi_actual).apply(rankdata, axis=0, method='average')
     
#Generate dictionary for all net loadings by variable and region
varContrib['USA']=zip(attr_names,attrib_contribution.tolist())

#############################################
##State Analysis     
#############################################
#Create New England conglomerate of states
US_All.loc[US_All.stateID.isin(['g23','g33','g25']), 'stateID'] = 'g23g33g25'

stateList = ['g23g33g25', 'g36','g51','g13','g17','g48','g29','g46','g06','g16']

for st in stateList:
    #Subset FEMA subregion
    stateData=US_All[US_All.stateID == st]

    # compute SoVI
    inputData = stateData.drop(['Geo_FIPS','stateID'], axis = 1, inplace = False)
    pca = SPSS_PCA(inputData, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=stateData.Geo_FIPS, columns=['sovi'])
    attrib_contribution = pca.weights_rot.sum(1)
    
    #sovi_alt_computation = (zinputs * attrib_contribution).sum(1) # this is just a check
    #sovi_alt_computation = pd.DataFrame(sovi_alt_computation, columns=['sovi'])
    #if not np.allclose(sovi_actual, sovi_alt_computation):
    #   raise Exception, "mismatch"
        
    ##Write attribute contribution output     
    #Generate dictionary for all net loadings by variable and region
    varContrib[st]=zip(attr_names,attrib_contribution.tolist())

#############################################
##Consolidate variable net contribs and ranks 
#############################################
netContribCols = varContrib.keys()

netContrib = pd.DataFrame(columns=netContribCols, index=attr_names)

for r in varContrib.keys():
    for name, value in varContrib[r]:
        netContrib.loc[name][r] = value

#variable rank using absolute value      
rankContrib = abs(netContrib).apply(rankdata, axis=0, method='average')
rankContrib = (28-rankContrib) + 1


combContrib = pd.DataFrame(columns=netContribCols, index=attr_names)
#can't think of a more elegant way to do this
for aRow in range(netContrib.shape[1]):
    for aCol in range(netContrib.shape[0]):
        combContrib.ix[aCol][aRow] = str(round(netContrib.ix[aCol][aRow], 2)) + ' (' + str(int(rankContrib.ix[aCol][aRow])) + ')'

#reorder table        
cols = ['USA', 'FEMA_1', 'g23g33g25', 
'FEMA_2', 'g36','FEMA_3', 'g51', 'FEMA_4', 'g13', 'FEMA_5', 'g17',
'FEMA_6', 'g48', 'FEMA_7', 'g29', 'FEMA_8', 'g46', 'FEMA_9', 'g06', 'FEMA_10', 
'g16']
combContrib = combContrib[cols]

#human readable variable names
desc = ['Median Age',
'Pop African-American (%)',
'Pop Native American (%)',
'Pop Asian (%)',
'Pop Hispanic (%)',
'Age Dependency (%)',
'Persons Per Housing Unit',
'Rental Housing (%)',
'Nursing Home Residents (%)',
'Pop Female (%)',
'Female-Headed Households (%)',
'Vacant Housing (%)',
'Per-Capita Income',
'English as Second Language (%)',
'Unemployment (%)',
'Poverty (%)',
'Mobile Homes (%)',
'Adults Completed <Grade 12 (%)',
'Female Employment (%)',
'Extractive Sector Employment (%)',
'Service Sector Employment (%)',
'Social Security Income (%)',
'No Automobile (%)',
'Children in Married Families (%)',
'Annual Income >$200K (%)',
'Median Rent',
'Median Home Value',
'Population Density']

#set descriptive names
combContrib.index = desc

#write out results
combContrib.to_csv(outPath+'/RegionVariableContribRank.csv')

#################
#Drop 1 Analysis
#################

geoLevels = US_All.Geo_FIPS
dropLevels = US_All.columns.drop(['Geo_FIPS', 'stateID'])
geoLabels = []
for _ in range(len(dropLevels)):
    geoLabels.extend(range(len(geoLevels)))
dropLabels = np.repeat(range(len(dropLevels)), len(geoLevels))

US_Drop1_Multi_Index = pd.MultiIndex(levels=[dropLevels, geoLevels], 
                                    labels=[dropLabels, geoLabels], 
                                    names=['DroppedVar', 'Geo_FIPS'])
                                    
US_Drop1_NetContrib = pd.DataFrame(index=dropLevels, columns=dropLevels)                     

US_SoVI_Drop1_Score = pd.DataFrame(index=US_Drop1_Multi_Index, columns=['sovi']) 


#Compute drop-one 
for j in dropLevels:
    US_dropj = US_All.drop([j,'Geo_FIPS', 'stateID'], axis = 1, inplace = False)
    pca = SPSS_PCA(US_dropj, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=geoLevels, columns=['sovi'])
    US_SoVI_Drop1_Score.loc[j, 'sovi'] = sovi_actual.values
    attrib_contribution = pd.DataFrame(data=pca.weights_rot.sum(1), index=US_dropj.columns)
    attrib_contribution = attrib_contribution.transpose()
    attrib_contribution.index = [j]
    US_Drop1_NetContrib.loc[j, attrib_contribution.columns] = attrib_contribution.values
    US_Drop1_NetContrib = US_Drop1_NetContrib.T #T so columns indexes dropped variable.

######################
##Compute spearman correlation
######################
#compute US base ranks

#rankContrib = (28-rankContrib) + 1
#compare Fema regions and US

#compare fema regions and states

#compare states and US

#compare drop 1 to full sovi

