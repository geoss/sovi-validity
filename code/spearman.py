import pandas as pd
from scipy.stats import spearmanr 

######################
##Compute spearman correlation
######################
##DO THIS USING THE USA AND STATE?FEMA REGION SPECIFIC SOVI

# compute US base ranks for all drop 1 correlations
# spearmanr(US_SoVI_Drop1_Score[])

# rankContrib = (28-rankContrib) + 1
# compare Fema regions and US

# compare fema regions and states

# compare states and US
# US_SoVI_Drop1_Score.
# compare drop 1 to full sovi

local_path = '/Users/becky/'
path = local_path + 'documents/sovi-validity/data/output/'

state_id = ['g51', 'g48', 'g36', 'g06', 'g13', 'g16', 'g17', 'g29', 'g46', 'g23g33g25']
rank = pd.read_csv(path + 'County_in_State_Rank.csv')
state = pd.read_csv(path + 'State_Sovi_Score.csv')
for ID in state_id:
    print ID
    st = state[state['state_id'] == ID]
    select = rank[rank['Geo_FIPS'].isin(st['Geo_FIPS'])]
    print "State: Region"
    print spearmanr(select['state_sovi_rank'], select['fema_region_sovi_rank'])
    print "State: US"
    print spearmanr(select['state_sovi_rank'], select['us_sovi_rank'])
    print "\n"

US_Sovi_Score = pd.read_csv(path + 'US_Sovi_Score.csv', index_col='Geo_FIPS')
FEMA_Region_Sovi_Score = pd.read_csv(path + 'FEMA_Region_Sovi_Score.csv', index_col='Geo_FIPS')
county_in_region_rank = pd.DataFrame(index=FEMA_Region_Sovi_Score.index,
                                    columns=['fema_region_sovi_rank', 'us_sovi_rank'])

regionList = ['FEMA_1','FEMA_2','FEMA_3','FEMA_4','FEMA_5','FEMA_6','FEMA_7','FEMA_8','FEMA_9','FEMA_10']

for region in regionList:
    x = FEMA_Region_Sovi_Score[FEMA_Region_Sovi_Score['fema_region'] == region]
    # get all counties in region and rank for us
    rg_cty_scores = US_Sovi_Score[US_Sovi_Score.index.isin(x.index)]
    county_in_region_rank.loc[rg_cty_scores.index, 'us_sovi_rank'] = abs(rg_cty_scores.sovi).rank(method='average', ascending=False)
    # get all counties in state and rank for fema region
    county_in_region_rank.loc[rg_cty_scores.index, 'fema_region_sovi_rank'] = abs(x.sovi).rank(method='average', ascending=False)

county_in_region_rank.to_csv(path + 'County_in_Region_Rank.csv')

for ID in regionList:
    print ID
    rg = FEMA_Region_Sovi_Score[FEMA_Region_Sovi_Score['fema_region'] == ID]
    rank = pd.read_csv(path + 'County_in_Region_Rank.csv')
    select = rank[rank['Geo_FIPS'].isin(rg.index)]
    print spearmanr(select['fema_region_sovi_rank'], select['us_sovi_rank'])
    
    
    