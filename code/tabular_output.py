#############################################
##Consolidate variable net contribs and ranks
#############################################
netContribCols = varContrib.keys()

netContrib = pd.DataFrame(columns=netContribCols, index=attr_names)

for r in varContrib.keys():
    for name, value in varContrib[r]:
        netContrib.loc[name][r] = value

netContrib = netContrib.reindex(netContrib.USA.abs().sort(inplace=False, ascending=False).index)

# reorder table
cols = ['USA', 'FEMA_1', 'g23g33g25',
        'FEMA_2', 'g36', 'FEMA_3', 'g51', 'FEMA_4', 'g13', 'FEMA_5', 'g17',
        'FEMA_6', 'g48', 'FEMA_7', 'g29', 'FEMA_8', 'g46', 'FEMA_9', 'g06', 'FEMA_10',
        'g16']
netContrib = netContrib[cols]

# variable rank using absolute value
rankContrib = abs(netContrib).apply(rankdata, axis=0, method='average')
rankContrib = (28 - rankContrib) + 1

combContrib = pd.DataFrame(columns=list(netContrib.columns), index=list(netContrib.index))
# can't think of a more elegant way to do this
for aRow in range(netContrib.shape[1]):
    for aCol in range(netContrib.shape[0]):
        combContrib.ix[aCol][aRow] = str(round(netContrib.ix[aCol][aRow], 2)) + ' (' + str(
            int(rankContrib.ix[aCol][aRow])) + ')'

# build list of varIDs and human readable names
# sort and use as index for conContrib
nameSort = [[name, hrname] for name, sign, sample, hrname in input_names]
nameSort = pd.DataFrame(nameSort)
nameSort.index = nameSort.loc[:, 0]
nameSort = nameSort.reindex(list(combContrib.index))

# set descriptive names
combContrib.index = list(nameSort.loc[:, 1])

# write out results
combContrib