import pandas as pd
from scipy.stats import spearmanr
import os
import sys

sys.path.insert(1, "./code")
from spss_pca import SPSS_PCA
import compute_sovis

pd.set_option("chained_assignment", None)

path = os.getcwd()
# path = os.path.dirname(os.getcwd()) # if running from the 'code' directory
outPath = os.path.join(path, 'data')
ipath = os.path.join(path, 'data', 'input')
spath = os.path.join(path, 'data', 'spatial')
opath = os.path.join(path, 'data', 'output')

state_id = ['g51', 'g48', 'g36', 'g06', 'g13', 'g16', 'g17', 'g29', 'g46', 'g23g33g25']
rank = compute_sovis.county_in_state_rank
state = compute_sovis.State_Sovi_Score
for ID in state_id:
    print ID
    st = state[state['state_id'] == ID]
    select = rank[rank.index.isin(st.index)]
    print "State: Region"
    print spearmanr(select['state_sovi_rank'], select['fema_region_sovi_rank'])
    print "State: US"
    print spearmanr(select['state_sovi_rank'], select['us_sovi_rank'])
    print "\n"

US_Sovi_Score = compute_sovis.US_Sovi_Score
FEMA_Region_Sovi_Score = compute_sovis.FEMA_Region_Sovi_Score
county_in_region_rank = pd.DataFrame(index=FEMA_Region_Sovi_Score.index,
                                     columns=['fema_region_sovi_rank', 'us_sovi_rank'])

regionList = ['FEMA_1', 'FEMA_2', 'FEMA_3', 'FEMA_4', 'FEMA_5', 'FEMA_6', 'FEMA_7', 'FEMA_8', 'FEMA_9', 'FEMA_10']

for region in regionList:
    x = FEMA_Region_Sovi_Score[FEMA_Region_Sovi_Score['fema_region'] == region]
    # get all counties in region and rank for us
    rg_cty_scores = US_Sovi_Score[US_Sovi_Score.index.isin(x.index)]
    county_in_region_rank.loc[rg_cty_scores.index, 'us_sovi_rank'] = abs(rg_cty_scores.sovi).rank(method='average',
                                                                                                  ascending=False)
    # get all counties in state and rank for fema region
    county_in_region_rank.loc[rg_cty_scores.index, 'fema_region_sovi_rank'] = abs(x.sovi).rank(method='average',
                                                                                               ascending=False)

county_in_region_rank.to_csv(path + 'County_in_Region_Rank.csv')

for ID in regionList:
    print ID
    rg = FEMA_Region_Sovi_Score[FEMA_Region_Sovi_Score['fema_region'] == ID]
    rank = pd.read_csv(path + 'County_in_Region_Rank.csv')
    select = rank[rank['Geo_FIPS'].isin(rg.index)]
    print spearmanr(select['fema_region_sovi_rank'], select['us_sovi_rank'])
