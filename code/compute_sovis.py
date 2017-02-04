# !/usr/bin/env

# This file calculates and exports five files
# US_Sovi_Score = a sovi analysis using the entire us outputs county score and rank
# FEMA_Region_Sovi_Score = a sovi analysis by fema region outputs county score and rank
# State_Sovi_Score = a sovi analysis by state for 10 states outputs county score and rank
# county_in_state_rank = a ranking of the counties of 10 states from the us, region level, and state level analysis
# variable_contributions = net contributions of each variable in each
# analysis above


import os
import sys
import pandas as pd
sys.path.insert(1, "./code")
from spss_pca import SPSS_PCA
import data_prep

pd.set_option("chained_assignment", None)

# Set paths to data
# local_path = '/Users/sspielman/'
# os.chdir(local_path + 'Dropbox/SoVI_var_wise_paper/code')
# path = local_path + '/Dropbox/SoVI_var_wise_paper'
# outPath = local_path + '/Dropbox/SoVI_var_wise_paper/data'
# ipath = local_path + "Dropbox/SoVI_var_wise_paper/data/input"
# spath = local_path + "Dropbox/SoVI_var_wise_paper/data/spatial"

path = os.getcwd()
# path = os.path.dirname(os.getcwd()) # if running from the 'code' directory
outPath = os.path.join(path, 'data')
ipath = os.path.join(path, 'data', 'input')
spath = os.path.join(path, 'data', 'spatial')

# copy db1 to new varname for clarity
US_All = data_prep.db1.copy()
US_All['Geo_FIPS'] = US_All.index.values
# US_All = pd.read_csv(path+'/data/input/sovi_inputs.csv')
# US_All.index = US_All.Geo_FIPS


# attribute name and expected influence on vulnerability
input_names = [['MEDAGE_ACS', 'pos', 'person', 'Median Age'],
               ['BLACK_ACS', 'pos', 'person', 'Pop African-American (%)'],
               ['QNATAM_ACS', 'pos', 'person', 'Pop Native American (%)'],
               ['QASIAN_ACS', 'pos', 'person', 'Pop Asian (%)'],
               ['QHISP_ACS', 'pos', 'person', 'Pop Hispanic (%)'],
               ['QAGEDEP_ACS', 'pos', 'person', 'Age Dependency (%)'],
               ['QPUNIT_ACS', 'pos', 'person', 'Persons Per Housing Unit'],
               ['PRENTER_ACS', 'pos', 'hu', 'Rental Housing (%)'],
               ['QNRRES_ACS', 'pos', 'person', 'Nursing Home Residents (%)'],
               ['QFEMALE_ACS', 'pos', 'person', 'Pop Female (%)'],
               ['QFHH_ACS', 'pos', 'hu', 'Female-Headed Households (%)'],
               ['QUNOCCHU_ACS', 'pos', 'hu', 'Vacant Housing (%)'],
               ['PERCAP_ALT', 'neg', 'person', 'Per-Capita Income'],
               ['QESL_ALT', 'pos', 'person', 'English as Second Language (%)'],
               ['QCVLUN', 'pos', 'person', 'Unemployment (%)'],
               ['QPOVTY', 'pos', 'person', 'Poverty (%)'],
               ['QMOHO', 'pos', 'hu', 'Mobile Homes (%)'],
               ['QED12LES_ALT', 'pos', 'person',
                   'Adults Completed <Grade 12 (%)'],
               ['QFEMLBR', 'pos', 'person', 'Female Employment (%)'],
               ['QEXTRCT_ALT', 'pos', 'person',
                   'Extractive Sector Employment (%)'],
               ['QSERV_ALT', 'pos', 'person', 'Service Sector Employment (%)'],
               ['QSSBEN', 'pos', 'hu', 'Social Security Income (%)'],
               ['QNOAUTO_ALT', 'pos', 'hu', 'No Automobile (%)'],
               ['QFAM', 'neg', 'person', 'Children in Married Families (%)'],
               ['QRICH200K', 'neg', 'hu', 'Annual Income >$200K (%)'],
               ['MDGRENT_ALT', 'neg', 'hu', 'Median Rent'],
               ['MHSEVAL_ALT', 'neg', 'hu', 'Median Home Value'],
               ['POPDENS', 'pos', 'person', 'Population Density']]

# Get attribute names
attr_names = [j[0] for j in input_names]
# cols = [c for c in US_All.columns if c.find('_SE') == -1]

attr_names.append('Geo_FIPS')
# US_All = US_All.dropna(axis=0) #two counties misisng data in state 15 and 48
US_All = US_All[attr_names]
US_All['stateID'] = US_All.Geo_FIPS.str.slice(0, 3, 1)
attr_names.remove('Geo_FIPS')

# ####Flipping Signs
# To ensure that each variable contributes as expected to the final Sovi
# Index following Eric Tate (2012?) we flip the signs of the input data.
for name, sign, sample, hrname in input_names:
    if sign == 'neg':
        US_All[name] = -US_All[name].values
    elif sign == 'pos':
        pass
    else:
        print("problem in flipping signs")
        raise

# Build FEMA subRegions Dict values= state ID's
FEMA_subs = dict()
FEMA_subs['FEMA_1'] = ['g23g33g25', 'g50', 'g09', 'g44']
FEMA_subs['FEMA_2'] = ['g36', 'g34']
FEMA_subs['FEMA_3'] = ['g42', 'g10', 'g11', 'g24', 'g51', 'g54']
FEMA_subs['FEMA_4'] = ['g21', 'g47', 'g37', 'g28', 'g01', 'g13', 'g45', 'g12']
FEMA_subs['FEMA_5'] = ['g27', 'g55', 'g26', 'g17', 'g18', 'g39']
FEMA_subs['FEMA_6'] = ['g35', 'g48', 'g40', 'g05', 'g22']
FEMA_subs['FEMA_7'] = ['g31', 'g19', 'g20', 'g29']
FEMA_subs['FEMA_8'] = ['g30', 'g38', 'g56', 'g46', 'g49', 'g08']
FEMA_subs['FEMA_9'] = ['g06', 'g32', 'g04']
FEMA_subs['FEMA_10'] = ['g53', 'g41', 'g16']

####################################
# DataFrames to hold US, fema region, and state level results
####################################

# Dict to hold variable loadings
# key will be [USA, Fema_region, stateid] depending on level of analysis
varContrib = {}

# National Score
US_Sovi_Score = pd.DataFrame(index=US_All.Geo_FIPS,
                             columns=['sovi', 'rank'])

# In the FEMA_Region_Sovi_Score data frame ranks are BY FEMA REGION.
# The data frame holds both the SOVI score and the county rank
# This means that there should be 10 counties with rank 1 (one for each
# FEMA Region)
FEMA_Region_Sovi_Score = pd.DataFrame(index=US_All.Geo_FIPS,
                                      columns=['sovi', 'rank', 'fema_region'])

# Create New England conglomerate of states
# These are the FIPS codes for the states wit hthe letter "g" appended
US_All.loc[US_All.stateID.isin(['g23', 'g33', 'g25']), 'stateID'] = 'g23g33g25'

# These are the states in the state level analysis
stateList = ['g23g33g25', 'g36', 'g51', 'g13',
             'g17', 'g48', 'g29', 'g46', 'g06', 'g16']

# In the State_Sovi_Score data frame ranks are BY STATE.
# The data frame holds both the SOVI score and the county rank
# This means that there should be 10 counties with rank 1 (one for each
# state in stateList)
State_Sovi_Score = pd.DataFrame(
    index=US_All.index[US_All.stateID.isin(stateList)],
    columns=['sovi', 'rank', 'state_id'])

#######################
# Compute National SoVI
#######################
# compute SoVI
inputData = US_All.drop(['Geo_FIPS', 'stateID'], axis=1, inplace=False)
pca = SPSS_PCA(inputData, reduce=True, varimax=True)
sovi_actual_us = pca.scores_rot.sum(1)
sovi_actual_us = pd.DataFrame(
    sovi_actual_us, index=US_All.Geo_FIPS, columns=['sovi'])
# rank
sovi_actual_us['rank'] = abs(sovi_actual_us).rank(
    method='average', ascending=False)
US_Sovi_Score.update(sovi_actual_us)

attrib_contribution_us = pca.weights_rot.sum(1)
# Generate dictionary for all net loadings by variable for US
varContrib['USA'] = zip(attr_names, attrib_contribution_us.tolist())


# quick check of ranks max should equal number of counties in US
try:
    US_Sovi_Score['rank'].max() == len(US_All)
except:
    print("error in ranking check")
    raise

# cleanup
del inputData
del sovi_actual_us
del attrib_contribution_us


######################
# FEMA REGION SOVI
######################
for i in FEMA_subs:

    # Subset FEMA subregion
    FEMARegionData = US_All[US_All['stateID'].isin(FEMA_subs[i])]

    # compute SoVI
    inputData = FEMARegionData.drop(
        ['Geo_FIPS', 'stateID'], axis=1, inplace=False)
    pca = SPSS_PCA(inputData, reduce=True, varimax=True)
    sovi_actual_fema = pca.scores_rot.sum(1)

    # load into df for merge
    sovi_actual_fema = pd.DataFrame(
        sovi_actual_fema, index=FEMARegionData.index, columns=['sovi'])
    # add fema region to df
    sovi_actual_fema['fema_region'] = i
    # rank
    sovi_actual_fema['rank'] = abs(sovi_actual_fema['sovi']).rank(
        method='average', ascending=False)

    FEMA_Region_Sovi_Score.update(sovi_actual_fema)

    attrib_contribution_fema = pca.weights_rot.sum(1)

    # Write attribute contribution output
    # Generate dictionary for all net loadings by variable and region
    varContrib[i] = zip(attr_names, attrib_contribution_fema.tolist())

# cleanup
del FEMARegionData
del inputData
del sovi_actual_fema
del attrib_contribution_fema

#############################################
# State Analysis
#############################################
for st in stateList:
    # Subset FEMA subregion
    stateData = US_All[US_All.stateID == st]

    # compute SoVI
    inputData = stateData.drop(['Geo_FIPS', 'stateID'], axis=1, inplace=False)
    pca = SPSS_PCA(inputData, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(
        sovi_actual, index=stateData.index, columns=['sovi'])
    sovi_actual['state_id'] = st
    # rank w/in state
    sovi_actual['rank'] = abs(sovi_actual['sovi']).rank(
        method='average', ascending=False)
    State_Sovi_Score.update(sovi_actual)
    attrib_contribution = pca.weights_rot.sum(1)
    varContrib[st] = zip(attr_names, attrib_contribution.tolist())

# cleanup
del stateData
del inputData
del sovi_actual
del attrib_contribution

###################################################
# Make Var Contributions Data Frame
###################################################
variable_contributions = pd.DataFrame(index=attr_names)
# for area in varContrib.iterkeys():
for area in varContrib.keys():
    variable_contributions[area] = [x for i, x in varContrib[area]]

##########################################################################
# For each county compute rank within state for US, state, and fema_region socis
##########################################################################

county_in_state_rank = pd.DataFrame(index=State_Sovi_Score.index,
                                    columns=['state_sovi_rank', 'fema_region_sovi_rank', 'us_sovi_rank'])

# for st in ['g23', 'g33', 'g25', 'g36', 'g51', 'g13', 'g17', 'g48', 'g29', 'g46', 'g06', 'g16']:
#     # get all counties in state and rank for us
#     st_cty_scores = US_Sovi_Score.loc[
#         [st in s for s in US_Sovi_Score.index], 'sovi']
#     county_in_state_rank.loc[st_cty_scores.index, 'us_sovi_rank'] = \
#         abs(st_cty_scores).rank(method='average', ascending=False)
#     # get all counties in state and rank for fema region
#     st_cty_scores = FEMA_Region_Sovi_Score.loc[
#         [st in s for s in FEMA_Region_Sovi_Score.index], 'sovi']
#     county_in_state_rank.loc[st_cty_scores.index, 'fema_region_sovi_rank'] = \
#         abs(st_cty_scores).rank(method='average', ascending=False)
#     # county rank in state only sovi
#     st_cty_scores = State_Sovi_Score.loc[
#         State_Sovi_Score['state_id'] == st, 'rank']
#     county_in_state_rank.loc[st_cty_scores.index,
#         'state_sovi_rank'] = st_cty_scores

for st in stateList:
    if st == 'g23g33g25':
        # get all counties in the three NE states and rank for us
        st_cty_scores = US_Sovi_Score.loc[
            ['g23' in s for s in US_Sovi_Score.index], 'sovi']
        st_cty_scores = st_cty_scores.append(
            US_Sovi_Score.loc[['g33' in s for s in US_Sovi_Score.index], 'sovi'])
        st_cty_scores = st_cty_scores.append(
            US_Sovi_Score.loc[['g25' in s for s in US_Sovi_Score.index], 'sovi'])

        county_in_state_rank.loc[st_cty_scores.index, 'us_sovi_rank'] = abs(st_cty_scores).rank(method='average', ascending=False)

        # get all counties in state and rank for fema region
        st_cty_scores = FEMA_Region_Sovi_Score.loc[
            ['g23' in s for s in FEMA_Region_Sovi_Score.index], 'sovi']
        st_cty_scores = st_cty_scores.append(FEMA_Region_Sovi_Score.loc[
                             ['g33' in s for s in FEMA_Region_Sovi_Score.index], 'sovi'])
        st_cty_scores = st_cty_scores.append(FEMA_Region_Sovi_Score.loc[
                             ['g25' in s for s in FEMA_Region_Sovi_Score.index], 'sovi'])

        county_in_state_rank.loc[st_cty_scores.index, 'fema_region_sovi_rank'] = abs(st_cty_scores).rank(method='average', ascending=False)

        st_cty_scores = State_Sovi_Score.loc[State_Sovi_Score['state_id'] == 'g23g33g25', 'rank']
        county_in_state_rank.loc[st_cty_scores.index, 'state_sovi_rank'] = st_cty_scores

    else:
        st_cty_scores = US_Sovi_Score.loc[[st in s for s in US_Sovi_Score.index], 'sovi']
        county_in_state_rank.loc[st_cty_scores.index, 'us_sovi_rank'] = abs(st_cty_scores).rank(method='average', ascending=False)
        # get all counties in state and rank for fema region
        st_cty_scores = FEMA_Region_Sovi_Score.loc[[st in s for s in FEMA_Region_Sovi_Score.index], 'sovi']
        county_in_state_rank.loc[st_cty_scores.index, 'fema_region_sovi_rank'] = abs(st_cty_scores).rank(method='average', ascending=False)
        # county rank in state only sovi
        st_cty_scores = State_Sovi_Score.loc[State_Sovi_Score['state_id'] == st, 'rank']
        county_in_state_rank.loc[st_cty_scores.index, 'state_sovi_rank'] = st_cty_scores


#####################################################
# OUTPUT TABLES
#####################################################
# US_Sovi_Score.to_csv(os.path.join(outPath,'output','US_Sovi_Score.csv')
# FEMA_Region_Sovi_Score.to_csv(os.path.join(outPath,'output','FEMA_Region_Sovi_Score.csv')
# State_Sovi_Score.to_csv(os.path.join(outPath,'output','State_Sovi_Score.csv')
# county_in_state_rank.to_csv(os.path.join(outPath,'output','County_in_State_Rank.csv')
# variable_contributions.to_csv(os.path.join(outPath,'output','variable_contributions.csv')

US_Sovi_Score.to_csv(os.path.join(outPath, 'output', 'US_Sovi_Score.csv'))
FEMA_Region_Sovi_Score.to_csv(os.path.join(
    outPath, 'output', 'FEMA_Region_Sovi_Score.csv'))
State_Sovi_Score.to_csv(os.path.join(
    outPath, 'output', 'State_Sovi_Score.csv'))
county_in_state_rank.to_csv(os.path.join(
    outPath, 'output', 'County_in_State_Rank.csv'))
variable_contributions.to_csv(os.path.join(
    outPath, 'output', 'variable_contributions.csv'))
