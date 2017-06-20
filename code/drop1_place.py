import os
import pandas as pd
# import geopandas as gpd
import pysal as ps
import numpy as np
from scipy.stats.mstats import zscore as ZSCORE
from scipy.stats import rankdata
from scipy.stats import spearmanr

# sovi compute script
# import sys
# sys.path.append(os.path.join(os.getcwd(),'code'))
from spss_pca import SPSS_PCA
# import compute_sovis

# outPath=compute_sovis.outPath

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

    # the data without no 1
    drop_no1=inputs.drop(drop)

    # preserve GEOIDs as an index
    # for computed SoVI
    geoLevels=drop_no1.Geo_FIPS

    #Compute drop "number one"
    pca = SPSS_PCA(drop_no1.drop(['Geo_FIPS', 'stateID'], axis = 1, inplace = False), reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=geoLevels, columns=['sovi'])
    drop1p_score = sovi_actual.values # SoVI score

    if netContrib is None: # no net contrib var's specified, compute county ranks

        # preserve the original county ranks
        # also for computing change in rank...
        orig_rank=scores.drop(drop)['rank']

        # add SoVI ranks for run
        drop1p_rank = pd.Series([i[0] for i in sovi_actual.values],index=geoLevels).rank(ascending=False)

        obs_rchg_drop1=pd.DataFrame({'orig_rank':orig_rank,'drop1p_rank':drop1p_rank},index=orig_rank.index)
        obs_rchg_drop1=obs_rchg_drop1.apply(lambda x: x.astype('int'),axis=1) # ensure all ints
        obs_rchg_drop1['rank_chg']=obs_rchg_drop1.orig_rank-obs_rchg_drop1.drop1p_rank

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

        nc_chg_drop1p=pd.DataFrame({'orig_rank':varRanks,'drop1p_rank':Drop1_NetContrib})
        nc_chg_drop1p=nc_chg_drop1p.apply(lambda x: x.astype('int'),axis=1) # ensure all ints
        nc_chg_drop1p['rank_chg']=nc_chg_drop1p.orig_rank-nc_chg_drop1p.drop1p_rank

        return nc_chg_drop1p

def rankChgTable(inputs,scores,obs_names,subset=None,top=5,cor=False,drop=1,verbose=True):

    dropany_result=dropAny(inputs=inputs,scores=scores,subset=subset,drop=drop)

    # ensure GEOID column in place
    # assumes missing GEOID values stored in df index
    if not 'geoFIPS' in dropany_result:
        dropany_result['geoFIPS']=dropany_result.index

    # merge dropany results with obs names
    rctab=dropany_result.merge(obs_names,on='geoFIPS')

    # print the spearman rank correlation if specified
    if cor:
        spearcor=spearmanr(rctab.drop1p_rank,rctab.orig_rank)
        if verbose:
            print("Spearman Rank Correlation: "+str(np.round(spearcor[0],5)),"\np-value: "+str(np.round(spearcor[1],4)))
            print('\n')

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
    d1rk=rctab[rctab.drop1p_rank<=top].ix[:,['geoFIPS','drop1p_rank','orig_rank','NAME']].sort_values('drop1p_rank')
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
        cor=spearmanr(drop_i.drop1p_rank,drop_i.orig_rank)
        cors.append(cor[0])

    return pd.Series(cors,index=geo_idx)

## function for plotting rank quantile moves

def rankQuantileMoves(inputs,scores,drop,subset=None,verbose=True):
    da=dropAny(inputs=inputs,scores=scores,subset=subset,drop=drop)
    if verbose:
        print(ps.Quantiles(da.orig_rank)) # quantile breaks key
        print('\n')
    r0=ps.Quantiles(da.orig_rank).yb
    r1=ps.Quantiles(da.drop1p_rank).yb
    moves_raw=pd.DataFrame({'r0':r0,'r1':r1}).groupby(['r0','r1']).size().unstack(fill_value=0)
    return np.round(moves_raw.apply(lambda x: x/sum(x),axis=1),2)
