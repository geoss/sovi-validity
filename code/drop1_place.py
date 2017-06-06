import os
import pandas as pd
import geopandas as gpd
import pysal as ps
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.mstats import zscore as ZSCORE
from scipy.stats import rankdata
from scipy.stats import spearmanr

# sovi compute script
# import sys
# sys.path.append(os.path.join(os.getcwd(),'code'))
from spss_pca import SPSS_PCA
import compute_sovis

outPath=compute_sovis.outPath

def dropAny(inputs,scores,drop=1,subset=None,netContrib=None,return_drop_rank=False):

    """
    inputs, the input variables, i.e. compute_sovis.USA_All

    scores, the SoVI outputs containing scores and ranks i.e. compute_sovis.USA_Sovi_Score

    subset: list of GEOIDs for subset (use for FEMA region or state)

    netContrib: i.e. compute_sovis.variable_ranks
    """

    if not subset is None:
        scores=scores[scores[scores.columns[len(scores.columns)-1]].astype('str').str.contains(subset)] # scores for region id
        inputs=inputs[inputs.index.isin(scores[scores.columns[len(scores.columns)-1]].index)] # inputs for subset

    # get the index of the
    # "most vulnerable" county
#     sovi_no1=scores[scores['rank']==scores['rank'].min()].index

#     if return_drop_rank:
#         return sovi_no1

    # the data without no 1
    drop_no1=inputs.drop(drop)

    # preserve GEOIDs as an index
    # for computed SoVI
    geoLevels=drop_no1.Geo_FIPS

    #Compute drop "number one"
    pca = SPSS_PCA(drop_no1.drop(['Geo_FIPS', 'stateID'], axis = 1, inplace = False), reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=geoLevels, columns=['sovi'])
    dropno1_score = sovi_actual.values # SoVI score

    if netContrib is None: # no net contrib var's specified, compute county ranks

        # preserve the original county ranks
        # also for computing change in rank...
        orig_rank=scores.drop(drop)['rank']

        # add SoVI ranks for run
        # any way to clean this up?
        dropno1_rank = pd.Series([i[0] for i in sovi_actual.values],index=geoLevels).rank(ascending=False)

        obs_rchg_drop1=pd.DataFrame({'orig_rank':orig_rank,'dropno1_rank':dropno1_rank},index=orig_rank.index)
        obs_rchg_drop1=obs_rchg_drop1.apply(lambda x: x.astype('int'),axis=1) # ensure all ints
        obs_rchg_drop1['rank_chg']=obs_rchg_drop1.orig_rank-obs_rchg_drop1.dropno1_rank

        return obs_rchg_drop1

    else: # net contribution

        #variable rank using absolute value
        rankContrib = abs(netContrib).apply(rankdata, axis=0, method='average')
        rankContrib = (28-rankContrib) + 1

        #Construct table to hold the results of the drop one analysis
        #Sort variable list based on importance rank.
        if not subset:
            varRanks = rankContrib['USA'].copy() #have to make a copy to sort index
            varRanks.sort('USA')
        else:
            varRanks = rankContrib[subset].copy() #have to make a copy to sort index
            varRanks.sort(subset)

        # recompute net contribution for drop no1
        Drop1_NetContrib = pd.Series(data=pca.weights_rot.sum(1), index=drop_no1.columns.drop(['Geo_FIPS', 'stateID']))
        Drop1_NetContrib = Drop1_NetContrib.transpose()
        Drop1_NetContrib=Drop1_NetContrib.convert_objects(convert_numeric=True)
        Drop1_NetContrib = Drop1_NetContrib.apply(lambda x: np.round(x, 2))
        Drop1_NetContrib = Drop1_NetContrib.rank(ascending=False)

        Drop1_NetContrib=Drop1_NetContrib[varRanks.index] # sort values by original index ranking

        nc_chg_dropno1=pd.DataFrame({'orig_rank':varRanks,'dropno1_rank':Drop1_NetContrib})
        nc_chg_dropno1=nc_chg_dropno1.apply(lambda x: x.astype('int'),axis=1) # ensure all ints
        nc_chg_dropno1['rank_chg']=nc_chg_dropno1.orig_rank-nc_chg_dropno1.dropno1_rank

        return nc_chg_dropno1

def rankChgTable(inputs,scores,obs_names,subset=None,top=5,cor=False,drop=1):

    dropany_result=dropAny(inputs=inputs,scores=scores,subset=subset,drop=drop)

    # ensure GEOID column in place
    # assumes missing GEOID values stored in df index
    if not 'geoFIPS' in dropany_result:
        dropany_result['geoFIPS']=dropany_result.index

    # merge dropany results with obs names
    rctab=dropany_result.merge(obs_names,on='geoFIPS')

    # print the spearman rank correlation if specified
    if cor:
        spearcor=spearmanr(rctab.dropno1_rank,rctab.orig_rank)
        print("Spearman Rank Correlation: "+str(np.round(spearcor[0],5)),"\np-value: "+str(np.round(spearcor[1],4)))

#     # get the dropped "number one" observation
#     if subset is None:
#         no1=obs_names[~obs_names.geoFIPS.isin(dropany_result.geoFIPS)]
#     else:
#         no1=obs_names[obs_names.geoFIPS==dropany(inputs=inputs,scores=scores,subset=subset,return_no1=True)[0]]
#     no1['orig_rank']=1

    # dropped obs rank
    drop_co=obs_names[obs_names.geoFIPS.str.contains(drop)]
    drop_co['orig_rank']=scores.ix[drop]['rank']

    # assemble table for original ranks
    orrk=rctab[rctab.orig_rank<=top].ix[:,['geoFIPS','orig_rank','NAME']]

    if int(drop_co.orig_rank)<=top: # append dropped obs to table if top-ranked
        orrk=drop_co.append(orrk)

    orrk=orrk.sort_values('orig_rank')
    orrk['Top_Orig']=orrk.NAME+" ("+orrk.orig_rank.astype('int').astype('str')+")"
    orrk.index=[i+1 for i in range(0,top)]
    orrk.ix[:,['NAME','Top_Orig']]

    # assemble table for dropany ranks
    d1rk=rctab[rctab.dropno1_rank<=top].ix[:,['geoFIPS','dropno1_rank','orig_rank','NAME']].sort_values('dropno1_rank')
    d1rk['Top_dropany']=d1rk.NAME+" ("+d1rk.orig_rank.astype('int').astype('str')+")"
    d1rk.index=[i+1 for i in range(0,top)]
    d1rk.ix[:,['NAME','Top_dropany']]

    # return the tables combined
    return pd.DataFrame({'All_Counties':orrk.Top_Orig,'Drop_1':d1rk.Top_dropany})

# wrap to a function
def dropCors(inputs,scores,subset=None):

    cors=[]

    if subset is None:
        geo_idx=scores.index.values
    else:
        geo_idx=scores[scores[scores.columns[len(scores.columns)-1]].astype('str').str.contains(subset)].index.values

    for i in geo_idx:
        drop_i=dropAny(inputs=inputs,scores=scores,subset=subset,drop=i)
        cor=spearmanr(drop_i.dropno1_rank,drop_i.orig_rank)
#         if cor[1]<0.1: # p-value threshold <0.1
#             cors.append(cor[0])
        cors.append(cor[0])

#     return cors
    return pd.Series(cors,index=geo_idx)

## function for plotting rank quantile moves

def rankQuantileMoves(inputs,scores,drop,subset=None):
    da=dropAny(inputs=inputs,scores=scores,subset=subset,drop=drop)
    print(ps.Quantiles(da.orig_rank)) # quantile breaks key
    r0=ps.Quantiles(da.orig_rank).yb
    r1=ps.Quantiles(da.dropno1_rank).yb
    moves_raw=pd.DataFrame({'r0':r0,'r1':r1}).groupby(['r0','r1']).size().unstack(fill_value=0)
    return np.round(moves_raw.apply(lambda x: x/sum(x),axis=1),2)

# df containing county names - no need for the geometries
# county_names=pd.DataFrame(gpd.read_file('data/spatial/USA_Counties_500k.shp').ix[:,['geoFIPS','NAME']],dtype='str')
county_names=pd.DataFrame(gpd.read_file('../data/spatial/USA_Counties_500k.shp').ix[:,['geoFIPS','NAME']],dtype='str') # loading from code folder


# ### States

# pd.unique(compute_sovis.State_Sovi_Score.state_id)


# ##### California

ca_cors=dropCors(compute_sovis.US_All,compute_sovis.State_Sovi_Score,'g06')
# ca_cors.describe()

cad=ca_cors[ca_cors==min(ca_cors)].index.values[0]
# county_names[county_names.geoFIPS.str.contains(cad)]

ca_rchg=rankChgTable(inputs=compute_sovis.US_All,scores=compute_sovis.State_Sovi_Score,obs_names=county_names,subset='g06',drop=cad,cor=True,top=10)
ca_rchg.to_csv(os.path.join(outPath,'output','ca_rank_change.csv'))

# rank quantile moves
ca_quint_moves=rankQuantileMoves(inputs=compute_sovis.US_All,scores=compute_sovis.State_Sovi_Score,subset='g06',drop=cad)
ca_quint_moves.to_csv(os.path.join(outPath,'output','ca_quint_moves.csv'))


# ##### FEMA 9: California and surrounding states (includes Hawaii)
f9_cors=dropCors(compute_sovis.US_All,compute_sovis.FEMA_Region_Sovi_Score,'FEMA_9')

# obs that decreases the correlation most when dropped
f9cd=f9_cors[f9_cors==min(f9_cors)].index.values[0]

f9_rchg=rankChgTable(inputs=compute_sovis.US_All,scores=compute_sovis.FEMA_Region_Sovi_Score,obs_names=county_names,subset='FEMA_9',drop=f9cd,cor=True,top=10)
f9_rchg.to_csv(os.path.join(outPath,'output','fema9_rank_change.csv'))

# rank quantile moves
f9_quint_moves=rankQuantileMoves(inputs=compute_sovis.US_All,scores=compute_sovis.FEMA_Region_Sovi_Score,subset='FEMA_9',drop=f9cd)
f9_quint_moves.to_csv(os.path.join(outPath,'output','fema9_quint_moves.csv'))


# ### Full USA
# NOT RUN - TIME INTENSIVE #
us_cors=dropCors(compute_sovis.US_All,compute_sovis.US_Sovi_Score)

# obs that decreases the correlation most when dropped
uscd=cors[cors==min(us_cors)].index.values[0]
county_names[county_names.geoFIPS.str.contains(uscd)]

us_rchg=rankChgTable(inputs=compute_sovis.US_All,scores=compute_sovis.US_Sovi_Score,obs_names=county_names,drop=uscd,cor=True,top=10)
us_rchg.to_csv(os.path.join(outPath,'output','usa_rank_change.csv'))

# rank quantile moves
us_rchg=rankQuantileMoves(inputs=compute_sovis.US_All,scores=compute_sovis.US_Sovi_Score,drop=uscd)
us_quint_moves.to_csv(os.path.join(outPath,'output','usa_quint_moves.csv'))
