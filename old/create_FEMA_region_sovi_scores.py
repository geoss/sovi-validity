#########################
##Script to create variable

import os as os
os.chdir('/Users/Seth/Dropbox/SoVI_var_wise_paper/code')

import numpy as np
import pandas as pd
from scipy.stats.mstats import zscore as ZSCORE
from scipy.stats import rankdata
from spss_pca import SPSS_PCA

path='/Users/Seth/Dropbox/SoVI_var_wise_paper'
outPath = '/Users/Seth/Dropbox/SoVI_var_wise_paper/data'

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
US_All = US_All.dropna(axis=0) #two counties misisng data in state 15 and 48
US_All = US_All[attr_names]
US_All['stateID'] = US_All.Geo_FIPS.str.slice(0,3,1)


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


for i in FEMA_subs:
    
    #Subset FEMA subregion
    FEMARegionData=US_All[US_All['stateID'].isin(FEMA_subs[i])]
    
    print "Subset",i+"...\n"

    # some data frames to organize the input data
    inputs = pd.DataFrame(index=FEMARegionData.Geo_FIPS)
    zinputs = pd.DataFrame(index=FEMARegionData.Geo_FIPS)

    
    # input data prep
    # --swap signs of the attributes expected to have a "negative" affect on vulnerability
    # --take z-score of data
    for name, sign, sample in input_names:
        if sign == 'neg':
            inputs[name] = -FEMARegionData[name].values
        elif sign == 'pos':
            inputs[name] = FEMARegionData[name].values
        else:
            raise Exception, "problem"
        zinputs[name] = ZSCORE(inputs[name].values)
    #zinputs_static = zinputs.copy(deep=True)
    
    # compute SoVI
    pca = SPSS_PCA(zinputs, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=FEMARegionData.Geo_FIPS, columns=['sovi'])
    
    # write SoVI to hard drive
    sovi_actual.to_csv(outPath+'/'+i+'_soviScores.csv')
    
    # share of each attribute that goes into SoVI
    # NOTE: since attributes are on different scales these shares cannot be
    #       compared across attributes
    attrib_contribution = pca.weights_rot.sum(1)
    sovi_alt_computation = (zinputs * attrib_contribution).sum(1) # this is just a check
    sovi_alt_computation = pd.DataFrame(sovi_alt_computation, columns=['sovi'])
    if not np.allclose(sovi_actual, sovi_alt_computation):
        raise Exception, "mismatch"
        
    ##Write attribute contribution output     
    #Generate dictionary for all net loadings by variable and region
    varContrib[i]=zip(attr_names,attrib_contribution.tolist())

#######################
##Compute National SoVI
#######################
zinputs = pd.DataFrame(index=US_All.Geo_FIPS)
inputs = pd.DataFrame(index=US_All.Geo_FIPS)
    
    # input data prep
    # --swap signs of the attributes expected to have a "negative" affect on vulnerability
    # --take z-score of data
for name, sign, sample in input_names:
    if sign == 'neg':
        inputs[name] = -US_All[name].values
    elif sign == 'pos':
        inputs[name] = US_All[name].values
    else:
        raise Exception, "problem"
    zinputs[name] = ZSCORE(US_All[name].values)

    
    # compute SoVI
pca = SPSS_PCA(zinputs, reduce=True, varimax=True)
sovi_actual = pca.scores_rot.sum(1)
sovi_actual = pd.DataFrame(sovi_actual, index=US_All.Geo_FIPS, columns=['sovi'])
    
    # write SoVI to hard drive
sovi_actual.to_csv(outPath+'/output/'+i+'_soviScores.csv')
    
    # share of each attribute that goes into SoVI
    # NOTE: since attributes are on different scales these shares cannot be
    #       compared across attributes
attrib_contribution = pca.weights_rot.sum(1)
sovi_alt_computation = (zinputs * attrib_contribution).sum(1) # this is just a check
sovi_alt_computation = pd.DataFrame(sovi_alt_computation, columns=['sovi'])
if not np.allclose(sovi_actual, sovi_alt_computation):
    raise Exception, "mismatch"
        
#Generate dictionary for all net loadings by variable and region
varContrib['USA']=zip(attr_names,attrib_contribution.tolist())

#############################################
##State Analysis     
#############################################
US_All.loc[US_All.stateID.isin(['g23','g33','g25']), 'stateID'] = 'g23g33g25'

stateList = ['g23g33g25', 'g36','g51','g13','g17','g48','g29','g46','g06','g16']

for st in stateList:
     #Subset FEMA subregion
    stateData=US_All[US_All.stateID == st]

    
    print "Subset",st+"...\n"

    # some data frames to organize the input data
    inputs = pd.DataFrame(index=stateData.Geo_FIPS)
    zinputs = pd.DataFrame(index=stateData.Geo_FIPS)

    
    # input data prep
    # --swap signs of the attributes expected to have a "negative" affect on vulnerability
    # --take z-score of data
    for name, sign, sample in input_names:
        if sign == 'neg':
            inputs[name] = -stateData[name].values
        elif sign == 'pos':
            inputs[name] = stateData[name].values
        else:
            raise Exception, "problem"
        zinputs[name] = ZSCORE(inputs[name].values)
    #zinputs_static = zinputs.copy(deep=True)
    
    # compute SoVI
    pca = SPSS_PCA(zinputs, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=stateData.Geo_FIPS, columns=['sovi'])
    
    # write SoVI to hard drive
    sovi_actual.to_csv(outPath+'/output/'+st+'_soviScores.csv')
    
    # share of each attribute that goes into SoVI
    # NOTE: since attributes are on different scales these shares cannot be
    #       compared across attributes
    attrib_contribution = pca.weights_rot.sum(1)
    sovi_alt_computation = (zinputs * attrib_contribution).sum(1) # this is just a check
    sovi_alt_computation = pd.DataFrame(sovi_alt_computation, columns=['sovi'])
    if not np.allclose(sovi_actual, sovi_alt_computation):
        raise Exception, "mismatch"
        
    ##Write attribute contribution output     
    #Generate dictionary for all net loadings by variable and region
    varContrib[st]=zip(attr_names,attrib_contribution.tolist())




#############################################
##Consolidate variable net contribs and ranks 
#############################################
netContribCols = [v for v in FEMA_subs.keys()] 
netContribCols.append('USA')
netContrib = pd.DataFrame(columns=netContribCols, index=attr_names)

for r in varContrib.keys():
    for name, value in varContrib[r]:
        netContrib.loc[name][r] = value

#variable rank using absolute value      
rankContrib = abs(netContrib).apply(rankdata, axis=0, method='average')
rankContrib = (28-rankContrib) + 1


combContrib = pd.DataFrame(columns=netContribCols, index=attr_names)
#cant htink of a more elegant way to do this
for aRow in range(netContrib.shape[1]):
    for aCol in range(netContrib.shape[0]):
        combContrib.ix[aCol][aRow] = str(round(netContrib.ix[aCol][aRow], 2)) + ' (' + str(int(rankContrib.ix[aCol][aRow])) + ')'

#reorder table        
#cols = ['USA', 'FEMA_1', 'FEMA_2', 'FEMA_3', 'FEMA_4', 'FEMA_5', 'FEMA_6', 'FEMA_7', 'FEMA_8', 'FEMA_9', 'FEMA_10']
#combContrib = combContrib[cols]

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

combContrib.to_csv(outPath+'/RegionVariableContribRank.csv')