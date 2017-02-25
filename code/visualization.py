<<<<<<< Updated upstream
=======
from bokeh.plotting import figures
import pandas as pd
from scipy.stats import spearmanr
import os
import sys

sys.path.insert(1, "./code")
from spss_pca import SPSS_PCA
import compute_sovis

# Construct table to hold the results of the drop one analysis
# Sort variable list based on importance rank.
USvarRanks = rankContrib.USA.copy()  # have to make a copy to sort index
USvarRanks.sort('USA')
dropLevels = USvarRanks.index

# build multindex
geoLevels = US_All.Geo_FIPS
geoLabels = []
for _ in range(len(dropLevels)):
    geoLabels.extend(range(len(geoLevels)))
dropLabels = np.repeat(range(len(dropLevels)), len(geoLevels))

US_Drop1_Multi_Index = pd.MultiIndex(levels=[dropLevels, geoLevels],
                                     labels=[dropLabels, geoLabels],
                                     names=['DroppedVar', 'Geo_FIPS'])

US_Drop1_NetContrib = pd.DataFrame(index=dropLevels, columns=dropLevels)

US_SoVI_Drop1_Score = pd.DataFrame(index=US_Drop1_Multi_Index, columns=['sovi'])

# Compute drop-one
for j in dropLevels:
    US_dropj = US_All.drop([j, 'Geo_FIPS', 'stateID'], axis=1, inplace=False)
    pca = SPSS_PCA(US_dropj, reduce=True, varimax=True)
    sovi_actual = pca.scores_rot.sum(1)
    sovi_actual = pd.DataFrame(sovi_actual, index=geoLevels, columns=['sovi'])
    US_SoVI_Drop1_Score.loc[j, 'sovi'] = sovi_actual.values
    attrib_contribution = pd.DataFrame(data=pca.weights_rot.sum(1), index=US_dropj.columns)
    # print(j +" " + str(np.isnan(attrib_contribution.values).sum()))
    attrib_contribution = attrib_contribution.transpose()
    attrib_contribution.index = [j]
    # print(attrib_contribution.loc[j,:])
    US_Drop1_NetContrib.loc[j, attrib_contribution.columns] = attrib_contribution.loc[j, :]  # .values

# Sort descriptive labels
USvarRanks = rankContrib.USA.copy()
USvarRanks.index = desc
USvarRanks.sort('USA')
US_Drop1_NetContrib.index = USvarRanks.index
US_Drop1_NetContrib.columns = USvarRanks.index

US_Drop1_NetContrib = US_Drop1_NetContrib.T  # T so columns indexes dropped variable.

# In[ ]:

US_Drop1_NetContrib = US_Drop1_NetContrib.convert_objects(convert_numeric=True)
US_Drop1_NetContrib = US_Drop1_NetContrib.apply(lambda x: np.round(x, 2))

# In[ ]:

get_ipython().magic(u'matplotlib inline')
sns.set_context("poster")

# Reorder and apply variable description to labels
# USvarRanks = rankContrib.USA.copy() #have to make a copy to sort index
# USvar
# dropLevels = USvarRanks.indexdesc

# plt.figure(figsize=(20, 16))
mask = np.isnan(US_Drop1_NetContrib)
sns.heatmap(US_Drop1_NetContrib, annot=True, linewidths=.25, vmin=-1, vmax=1, annot_kws={"size": 7})
>>>>>>> Stashed changes
