# !/usr/bin/env

# This script loads and prepares Census and ACS data.
# Outputs a CSV file that can be used for the construction of SoVI
# Interpolates missing values using nearby values
# Calculates SE for compositve variables

import os
import pandas as pd
import pysal as ps
import numpy as np

pd.set_option("chained_assignment", None)

# # Set paths to data
# local_path = '/Users/sspielman/'
# os.chdir(local_path + 'Dropbox/SoVI_var_wise_paper/code')
# path = local_path + '/Dropbox/SoVI_var_wise_paper'
# outPath = local_path + '/Dropbox/SoVI_var_wise_paper/data'
# ipath = local_path + "Dropbox/SoVI_var_wise_paper/data/input"
# spath = local_path + "Dropbox/SoVI_var_wise_paper/data/spatial"

path = os.getcwd()
# path = os.path.dirname(os.getcwd()) # if running from the 'code' directory
outPath=os.path.join(path,'data')
ipath = os.path.join(path,'data','input')
spath = os.path.join(path,'data','spatial')
#
# functions fot the calculation of SE


def se_sum(*ses):
    df_temp = pd.DataFrame(list(ses))
    df_temp = df_temp.T
    df_temp = np.square(df_temp)
    df_temp = df_temp.sum(1)
    return np.sqrt(df_temp)


# SE of a ratio


def se_ratio(est, estd, sen, sed):
    sen2 = np.square(sen)
    sed2 = np.square(sed)
    est2 = np.square(est)
    num = np.sqrt(sen2 + (est2 * sed2))
    return num / estd


# SE of a proprotion


def se_prop(est, estd, sen, sed):
    sen2 = np.square(sen)
    sed2 = np.square(sed)
    est2 = np.square(est)
    num = sen2 - (est2 * sed2)
    num_alt = sen2 + (est2 * sed2)
    problems = num <= 0
    num[problems] = num_alt[problems]
    num = np.sqrt(num)
    return num / estd


# unit test for equivalency between original and constructed variables


def equal_test(orig, alt):
    if np.equal(orig, alt).sum() != db.shape[0]:
        if (db.shape[0] - np.equal(orig, alt).sum()) \
                == np.isnan(orig).sum() == np.isnan(alt).sum():
            pass
        else:
            print("problem in equal test")
            raise


# define data types
make_strings = {'Geo_FIPS': object, 'Geo_STATE': object, 'Geo_COUNTY': object,
                'Geo_TRACT': object, 'Geo_CBSA': object, 'Geo_CSA': object}

# load data
# JOE: I had to change the encoding to latin-1 to avoid hitting a UTF-8 error
acs = pd.read_csv(os.path.join(ipath, 'sovi_acs.csv'),
                  dtype=make_strings, skiprows=1,encoding='latin-1')
census = pd.read_csv(os.path.join(ipath, 'sovi_decennial.csv'),
                     dtype=make_strings, skiprows=1,encoding='latin-1')
acs_samp = pd.read_csv(os.path.join(ipath, 'sovi_acs_sampSize.csv'),
                       dtype=make_strings, skiprows=1,encoding='latin-1')

# format FIPS
acs.index = 'g' + acs.Geo_FIPS
census.index = 'g' + census.Geo_FIPS
acs_samp.index = 'g' + acs_samp.Geo_FIPS

# merge files
db = census
db = db.join(acs, rsuffix='_acs')
db = db.join(acs_samp, rsuffix='_acsSamp')

# if available add supplmentary data
try:
    census_sup1 = pd.read_csv(os.path.join(ipath, 'sovi_decennial_sup1.csv'),
        dtype=make_strings,skiprows=1,encoding='latin-1')
    census_sup1.index = 'g' + census_sup1.Geo_FIPS
    db = db.join(census_sup1, rsuffix='_decSup1')
except:
    print("no supplementary decennial data")
try:
    acs_sup1 = pd.read_csv(os.path.join(spath, 'sovi_acs_sup1.csv'),
        dtype=make_strings,skiprows=1,encoding='latin-1')
    acs_sup1.index = 'g' + acs_sup1.Geo_FIPS
    db = db.join(acs_sup1, rsuffix='_acsSup1')
except:
    print("did not pull supplementary ACS data - A")
try:
    acs_sup2 = pd.read_csv(os.path.join(ipath, 'sovi_acs_kids.csv'),
                           dtype=make_strings, skiprows=1,encoding='latin-1')
    acs_sup2.index = 'g' + acs_sup2.Geo_FIPS
    db = db.join(acs_sup2, rsuffix='_acsSup2')
except:
    print("did not pull supplementary ACS data - B")

# drop Puerto Rico (sorry PR!)
db = db[db.Geo_STATE != '72']

# define SE columns
se_cols = [i for i in db.columns if i[-1] == 's' and i[0] == 'A']
db[se_cols] *= (1.65 / 1.645)

# calculate weights matrix
w = ps.queen_from_shapefile(os.path.join(spath, 'USA_Counties_500k.shp'),
                            idVariable='geoFIPS')
w.transform = 'R'

# output dataframe
db1 = pd.DataFrame(index=db.index)

# Decennial variables (original)
db1['MEDAGE'] = db.SF1_P0130001
db1['BLACK'] = (db.SF1_P0030003 * 1.) / db.SF1_P0010001
db1['QNATAM'] = (db.SF1_P0030004 * 1.) / db.SF1_P0010001
db1['QASIAN'] = (db.SF1_P0030005 * 1.) / db.SF1_P0010001
db1['QHISP'] = (db.SF1_P0040003 * 1.) / db.SF1_P0010001
db1['QAGEDEP'] = ((db.SF1_P0120003 + db.SF1_P0120027 + db.SF1_P0120020 +
                   db.SF1_P0120021 + db.SF1_P0120022 + db.SF1_P0120023 +
                   db.SF1_P0120024 + db.SF1_P0120025 + db.SF1_P0120044 +
                   db.SF1_P0120045 + db.SF1_P0120046 + db.SF1_P0120047 +
                   db.SF1_P0120048 + db.SF1_P0120049) * 1.) / db.SF1_P0010001
db1['PPUNIT'] = db.SF1_H0100001 / (db.SF1_H0030002 * 1.)
db1['PRENTER'] = (db.SF1_H0040004 * 1.) / db.SF1_H0010001
db1['QNRRES'] = (db.SF1_P0420005 * 1.) / db.SF1_P0010001
db1['QFEMALE'] = (db.SF1_P0120026 * 1.) / db.SF1_P0010001
db1['QFHH'] = (db.SF1_P0190014 * 1.) / db.SF1_P0180001
db1['QUNOCCHU'] = ((db.SF1_H0010001 - db.SF1_H0030002) * 1.) / db.SF1_H0010001

# Decennial variables (alternatives)
db1['BLACK_ALT'] = (db.SF1_P0050004 * 1.) / db.SF1_P0010001  # exclude hispanic
db1['QNATAM_ALT'] = (db.SF1_P0050005 * 1.) / \
                    db.SF1_P0010001  # exclude hispanic
db1['QASIAN_ALT'] = (db.SF1_P0050006 * 1.) / \
                    db.SF1_P0010001  # exclude hispanic
db1['QNRRES_ALT'] = (db.SF1_P0430023 + db.SF1_P0430054 * 1.) / \
                    db.SF1_P0010001  # 65 and over living in group quarters
# same value, simplified computation
db1['QUNOCCHU_ALT'] = (db.SF1_H0030003 * 1.) / db.SF1_H0030001

# Decennial variables (using ACS data and alternative formulations)
db1['MEDAGE_ACS'] = db.ACS12_5yr_B01002001
db1['BLACK_ACS'] = db.ACS12_5yr_B03002004 / (db.ACS12_5yr_B03002001 * 1.)
db1['QNATAM_ACS'] = db.ACS12_5yr_B03002005 / (db.ACS12_5yr_B03002001 * 1.)
db1['QASIAN_ACS'] = db.ACS12_5yr_B03002006 / (db.ACS12_5yr_B03002001 * 1.)
db1['QHISP_ACS'] = db.ACS12_5yr_B03002012 / (db.ACS12_5yr_B03002001 * 1.)
db1['QAGEDEP_ACS'] = (db.ACS12_5yr_B06001002 +
                      db.ACS12_5yr_B09020001) / (db.ACS12_5yr_B01003001 * 1.)
db1['QPUNIT_ACS'] = db.ACS12_5yr_B25008001 / (db.ACS12_5yr_B25002002 * 1.)
db1['PRENTER_ACS'] = db.ACS12_5yr_B25003003 / (db.ACS12_5yr_B25002001 * 1.)
db1['QNRRES_ACS'] = db.ACS12_5yr_B09020021 / (db.ACS12_5yr_B01003001 * 1.)
db1['QFEMALE_ACS'] = db.ACS12_5yr_B01001026 / (db.ACS12_5yr_B01003001 * 1.)
db1['QFHH_ACS'] = db.ACS12_5yr_B11001006 / (db.ACS12_5yr_B11001001 * 1.)
db1['QUNOCCHU_ACS'] = db.ACS12_5yr_B25002003 / (db.ACS12_5yr_B25002001 * 1.)

# ACS variables (original)
db1['PERCAP'] = db.ACS12_5yr_B19025001 / (db.ACS12_5yr_B01003001 * 1.)
db1['QESL'] = ((db.ACS12_5yr_B16004029 + db.ACS12_5yr_B16004030 +
                db.ACS12_5yr_B16004034 + db.ACS12_5yr_B16004035 +
                db.ACS12_5yr_B16004039 + db.ACS12_5yr_B16004040 +
                db.ACS12_5yr_B16004044 + db.ACS12_5yr_B16004045 +
                db.ACS12_5yr_B16004051 + db.ACS12_5yr_B16004052 +
                db.ACS12_5yr_B16004056 + db.ACS12_5yr_B16004057 +
                db.ACS12_5yr_B16004061 + db.ACS12_5yr_B16004062 +
                db.ACS12_5yr_B16004066 + db.ACS12_5yr_B16004067) * 1.) / \
              ((db.ACS12_5yr_B16004024 + db.ACS12_5yr_B16004046) -
               (db.ACS12_5yr_B16004025 + db.ACS12_5yr_B16004047))
db1.QESL = db1.QESL.replace([np.inf, -np.inf, np.nan], 0)
db1.QESL = db1.QESL.replace([np.inf, -np.inf], 0)
db1['QCVLUN'] = ((db.ACS12_5yr_B23022025 + db.ACS12_5yr_B23022049) * 1.) / \
                db.ACS12_5yr_B23022001
db1['QPOVTY'] = (db.ACS12_5yr_B17021002 * 1.) / db.ACS12_5yr_B17021001
db1['QMOHO'] = (db.ACS12_5yr_B25024010 * 1.) / db.ACS12_5yr_B25024001
db1['QED12LES'] = ((db.ACS12_5yr_B15002003 + db.ACS12_5yr_B15002004 +
                    db.ACS12_5yr_B15002005 + db.ACS12_5yr_B15002006 +
                    db.ACS12_5yr_B15002007 + db.ACS12_5yr_B15002008 +
                    db.ACS12_5yr_B15002009 + db.ACS12_5yr_B15002010 +
                    db.ACS12_5yr_B15002020 + db.ACS12_5yr_B15002021 +
                    db.ACS12_5yr_B15002022 + db.ACS12_5yr_B15002023 +
                    db.ACS12_5yr_B15002024 + db.ACS12_5yr_B15002025 +
                    db.ACS12_5yr_B15002026 + db.ACS12_5yr_B15002027) * 1.) / \
                  db.ACS12_5yr_B15002001
db1['QFEMLBR'] = (db.ACS12_5yr_C24010038 * 1.) / db.ACS12_5yr_C24010001
db1['QEXTRCT'] = ((db.ACS12_5yr_C24030003 + db.ACS12_5yr_C24030030) * 1.) / \
                 db.ACS12_5yr_C24030001
db1['QSERV'] = ((db.ACS12_5yr_C24010019 + db.ACS12_5yr_C24010055) * 1.) / \
               db.ACS12_5yr_C24010001
db1['QSSBEN'] = (db.ACS12_5yr_B19055002 * 1.) / db.ACS12_5yr_B19055001
db1['QNOAUTO'] = ((db.ACS12_5yr_B25044003 + db.ACS12_5yr_B25044010) * 1.) / \
                 db.ACS12_5yr_B25044001
db1['QFAM'] = (db.ACS12_5yr_B09002002 * 1.) / db.ACS12_5yr_B09002001
db1.QFAM = db1.QFAM.replace([np.inf, -np.inf, np.nan], 0)
db1['QRICH200K'] = (db.ACS12_5yr_B19001017 * 1.) / db.ACS12_5yr_B11001001

# ACS variables (alternatives)
# HH income divided by persons in HHs
db1['PERCAP_ALT'] = db.ACS12_5yr_B19025001 / (db.ACS12_5yr_B25008001 * 1.)

# 5 and older who don't speak English very well
db1['QESL_ALT'] = ((db.ACS12_5yr_B06007005 + db.ACS12_5yr_B06007008) * 1.) / \
                  db.ACS12_5yr_B06007001

# same value, simplified computation
db1['QED12LES_ALT'] = (db.ACS12_5yr_B16010002 * 1.) / db.ACS12_5yr_B16010001

# same value, simplified computation
db1['QEXTRCT_ALT'] = (db.ACS12_5yr_C24050002 * 1.) / db.ACS12_5yr_C24050001

# same value, simplified computation
db1['QSERV_ALT'] = (db.ACS12_5yr_C24050029 * 1.) / db.ACS12_5yr_C24050001

# same value, simplified computation
db1['QNOAUTO_ALT'] = (db.ACS12_5yr_B08201002 * 1.) / db.ACS12_5yr_B08201001

# the original computed the median by hand so is not included
db1['MDGRENT_ALT'] = db.ACS12_5yr_B25064001

# the original computed the median by hand so is not included
db1['MHSEVAL_ALT'] = db.ACS12_5yr_B25077001

# I didn't understand QURBRURX
db1['POPDENS'] = db.ACS12_5yr_B01003001 / (db.SE_T02A_002 * 1.)

# if no home value, assign the spatial lag of the estimate and SE
homeval = db1['MHSEVAL_ALT'].copy()
homeval_se = db.ACS12_5yr_B25077001s.copy()
dbf = ps.open(os.path.join(spath, 'USA_Counties_500k.dbf'))

# Rename dbf GEOIDs to match homeval
geoid = dbf.by_col('geoFIPS')

shp_fips = pd.DataFrame(dbf.by_col('geoFIPS'), index=geoid)
shp_fips = shp_fips.join(homeval)
shp_fips = shp_fips.join(homeval_se)
shp_fips['MHSEVAL_ALT_LAG'] = ps.lag_spatial(w, shp_fips.MHSEVAL_ALT)
shp_fips['MHSEVAL_ALT_LAG_SE'] = ps.lag_spatial(w, shp_fips.ACS12_5yr_B25077001s)

mh = shp_fips.ix[shp_fips.MHSEVAL_ALT_LAG == 0].MHSEVAL_ALT.tolist()

# Reassign values to MHSEVAL_ALT_LAG
shp_fips.ix[shp_fips.MHSEVAL_ALT_LAG == 0, 'MHSEVAL_ALT_LAG'] = mh

# Reassign missing standard error values
mhs = shp_fips.ix[shp_fips.MHSEVAL_ALT_LAG_SE == 0].ACS12_5yr_B25077001s.tolist()
shp_fips.ix[shp_fips.MHSEVAL_ALT_LAG_SE == 0, 'MHSEVAL_ALT_LAG_SE'] = mhs

# Get rid of nan values - reassign MHSEVAL_ALT(_SE)
shp_fips.MHSEVAL_ALT_LAG[np.isnan(shp_fips.MHSEVAL_ALT_LAG)] = \
    shp_fips.MHSEVAL_ALT[np.isnan(shp_fips.MHSEVAL_ALT_LAG)]  # replace NA with lag
shp_fips.MHSEVAL_ALT_LAG_SE[np.isnan(shp_fips.MHSEVAL_ALT_LAG_SE)] = \
    shp_fips.ACS12_5yr_B25077001s[np.isnan(shp_fips.MHSEVAL_ALT_LAG_SE)]  # replace NA with lag

db1['MHSEVAL_ALT_LAG'] = shp_fips['MHSEVAL_ALT_LAG']
db1['MHSEVAL_ALT_LAG_SE'] = shp_fips['MHSEVAL_ALT_LAG_SE']
db1.MHSEVAL_ALT[np.isnan(db1.MHSEVAL_ALT)] = db1.MHSEVAL_ALT_LAG[np.isnan(db1.MHSEVAL_ALT)]
# note: the lagged SE is pushed to the final column in the SE section below

#############################

# Decennial standard errors (using ACS data and alternative formulations)
db1['MEDAGE_ACS_SE'] = db.ACS12_5yr_B01002001s

db1['BLACK_ACS_SE'] = se_prop(db1.BLACK_ACS, db.ACS12_5yr_B03002001,
                              db.ACS12_5yr_B03002004s, db.ACS12_5yr_B03002001s)
db1['QNATAM_ACS_SE'] = se_prop(db1.QNATAM_ACS, db.ACS12_5yr_B03002001,
                               db.ACS12_5yr_B03002005s, db.ACS12_5yr_B03002001s)
db1['QASIAN_ACS_SE'] = se_prop(db1.QASIAN_ACS, db.ACS12_5yr_B03002001,
                               db.ACS12_5yr_B03002006s, db.ACS12_5yr_B03002001s)
db1['QHISP_ACS_SE'] = se_prop(db1.QHISP_ACS, db.ACS12_5yr_B03002001,
                              db.ACS12_5yr_B03002012s, db.ACS12_5yr_B03002001s)

QAGEDEP_ACS_sen = se_sum(db.ACS12_5yr_B06001002s, db.ACS12_5yr_B09020001s)
db1['QAGEDEP_ACS_SE'] = se_prop(db1.QAGEDEP_ACS, db.ACS12_5yr_B01003001,
                                QAGEDEP_ACS_sen, db.ACS12_5yr_B01003001s)

db1['QPUNIT_ACS_SE'] = se_ratio(db1.QPUNIT_ACS, db.ACS12_5yr_B25002002,
                                db.ACS12_5yr_B25008001s, db.ACS12_5yr_B25002002s)
db1['PRENTER_ACS_SE'] = se_prop(db1.PRENTER_ACS, db.ACS12_5yr_B25002001,
                                db.ACS12_5yr_B25003003s, db.ACS12_5yr_B25002001s)
db1['QNRRES_ACS_SE'] = se_prop(db1.QNRRES_ACS, db.ACS12_5yr_B01003001,
                               db.ACS12_5yr_B09020021s, db.ACS12_5yr_B01003001s)
db1['QFEMALE_ACS_SE'] = se_prop(db1.QFEMALE_ACS, db.ACS12_5yr_B01003001,
                                db.ACS12_5yr_B01001026s, db.ACS12_5yr_B01003001s)
db1['QFHH_ACS_SE'] = se_prop(db1.QFHH_ACS, db.ACS12_5yr_B11001001,
                             db.ACS12_5yr_B11001006s, db.ACS12_5yr_B11001001s)
db1['QUNOCCHU_ACS_SE'] = se_prop(db1.QUNOCCHU_ACS, db.ACS12_5yr_B25002001,
                                 db.ACS12_5yr_B25002003s, db.ACS12_5yr_B25002001s)

# ACS standard errors (original)
db1['PERCAP_SE'] = se_ratio(db1.PERCAP, db.ACS12_5yr_B01003001,
                            db.ACS12_5yr_B19025001s, db.ACS12_5yr_B01003001s)

QESL_sen = se_sum(db.ACS12_5yr_B16004029s, db.ACS12_5yr_B16004030s,
                  db.ACS12_5yr_B16004034s, db.ACS12_5yr_B16004035s,
                  db.ACS12_5yr_B16004039s, db.ACS12_5yr_B16004040s,
                  db.ACS12_5yr_B16004044s, db.ACS12_5yr_B16004045s,
                  db.ACS12_5yr_B16004051s, db.ACS12_5yr_B16004052s,
                  db.ACS12_5yr_B16004056s, db.ACS12_5yr_B16004057s,
                  db.ACS12_5yr_B16004061s, db.ACS12_5yr_B16004062s,
                  db.ACS12_5yr_B16004066s, db.ACS12_5yr_B16004067s)
QESL_sed = se_sum(db.ACS12_5yr_B16004024s, db.ACS12_5yr_B16004046s,
                  db.ACS12_5yr_B16004025s, db.ACS12_5yr_B16004047s)
db1['QESL_SE'] = se_prop(db1.QESL, (db.ACS12_5yr_B16004024 +
                                    db.ACS12_5yr_B16004046) -
                         (db.ACS12_5yr_B16004025 + db.ACS12_5yr_B16004047),
                         QESL_sen, QESL_sed)
db1.QESL_SE = db1.QESL_SE.replace([np.inf, -np.inf], 0)
db1.QESL_SE[db1.QESL == 0] = 0

QCVLUN_sen = se_sum(db.ACS12_5yr_B23022025s, db.ACS12_5yr_B23022049s)
db1['QCVLUN_SE'] = se_prop(db1.QCVLUN, db.ACS12_5yr_B23022001,
                           QCVLUN_sen, db.ACS12_5yr_B23022001s)

db1['QPOVTY_SE'] = se_prop(db1.QPOVTY, db.ACS12_5yr_B17021001,
                           db.ACS12_5yr_B17021002s, db.ACS12_5yr_B17021001s)
db1['QMOHO_SE'] = se_prop(db1.QMOHO, db.ACS12_5yr_B25024001,
                          db.ACS12_5yr_B25024010s, db.ACS12_5yr_B25024001s)

QED12LES_sen = se_sum(db.ACS12_5yr_B15002003s, db.ACS12_5yr_B15002004s,
                      db.ACS12_5yr_B15002005s, db.ACS12_5yr_B15002006s,
                      db.ACS12_5yr_B15002007s, db.ACS12_5yr_B15002008s,
                      db.ACS12_5yr_B15002009s, db.ACS12_5yr_B15002010s,
                      db.ACS12_5yr_B15002020s, db.ACS12_5yr_B15002021s,
                      db.ACS12_5yr_B15002022s, db.ACS12_5yr_B15002023s,
                      db.ACS12_5yr_B15002024s, db.ACS12_5yr_B15002025s,
                      db.ACS12_5yr_B15002026s, db.ACS12_5yr_B15002027s)
db1['QED12LES_SE'] = se_prop(db1.QED12LES, db.ACS12_5yr_B15002001,
                             QED12LES_sen, db.ACS12_5yr_B15002001s)

db1['QFEMLBR_SE'] = se_prop(db1.QFEMLBR, db.ACS12_5yr_C24010001,
                            db.ACS12_5yr_C24010038s, db.ACS12_5yr_C24010001s)

QEXTRCT_sen = se_sum(db.ACS12_5yr_C24030003s, db.ACS12_5yr_C24030030s)
db1['QEXTRCT_SE'] = se_prop(db1.QEXTRCT, db.ACS12_5yr_C24030001,
                            QEXTRCT_sen, db.ACS12_5yr_C24030001s)

QSERV_sen = se_sum(db.ACS12_5yr_C24010019s, db.ACS12_5yr_C24010055s)
db1['QSERV_SE'] = se_prop(db1.QSERV, db.ACS12_5yr_C24010001,
                          QSERV_sen, db.ACS12_5yr_C24010001s)

db1['QSSBEN_SE'] = se_prop(db1.QSSBEN, db.ACS12_5yr_B19055001,
                           db.ACS12_5yr_B19055002s, db.ACS12_5yr_B19055001s)

QNOAUTO_sen = se_sum(db.ACS12_5yr_B25044003s, db.ACS12_5yr_B25044010s)
db1['QNOAUTO_SE'] = se_prop(db1.QNOAUTO, db.ACS12_5yr_B25044001,
                            QNOAUTO_sen, db.ACS12_5yr_B25044001s)

db1['QFAM_SE'] = se_prop(db1.QFAM, db.ACS12_5yr_B09002001,
                         db.ACS12_5yr_B09002002s, db.ACS12_5yr_B09002001s)
db1.QFAM_SE = db1.QFAM_SE.replace([np.inf, -np.inf], 0)

db1['QRICH200K_SE'] = se_prop(db1.QRICH200K, db.ACS12_5yr_B11001001,
                              db.ACS12_5yr_B19001017s, db.ACS12_5yr_B11001001s)

#############################

# ACS standard errors (alternatives)
db1['PERCAP_ALT_SE'] = se_ratio(db1.PERCAP_ALT, db.ACS12_5yr_B25008001,
                                db.ACS12_5yr_B19025001s,
                                db.ACS12_5yr_B25008001s)

QESL_ALT_sen = se_sum(db.ACS12_5yr_B06007005s, db.ACS12_5yr_B06007008s)
db1['QESL_ALT_SE'] = se_prop(db1.QESL_ALT, db.ACS12_5yr_B06007001,
                             QESL_ALT_sen, db.ACS12_5yr_B06007001s)

db1['QED12LES_ALT_SE'] = se_prop(db1.QED12LES_ALT, db.ACS12_5yr_B16010001,
                                 db.ACS12_5yr_B16010002s, db.ACS12_5yr_B16010001s)
db1['QEXTRCT_ALT_SE'] = se_prop(db1.QEXTRCT_ALT, db.ACS12_5yr_C24050001,
                                db.ACS12_5yr_C24050002s, db.ACS12_5yr_C24050001s)
db1['QSERV_ALT_SE'] = se_prop(db1.QSERV_ALT, db.ACS12_5yr_C24050001,
                              db.ACS12_5yr_C24050029s, db.ACS12_5yr_C24050001s)
db1['QNOAUTO_ALT_SE'] = se_prop(db1.QNOAUTO_ALT, db.ACS12_5yr_B08201001,
                                db.ACS12_5yr_B08201002s, db.ACS12_5yr_B08201001s)
db1['MDGRENT_ALT_SE'] = db.ACS12_5yr_B25064001s
db1['MHSEVAL_ALT_SE'] = db.ACS12_5yr_B25077001s
db1.MHSEVAL_ALT_SE[np.isnan(db1.MHSEVAL_ALT)] = db1.MHSEVAL_ALT_LAG_SE[
    np.isnan(db1.MHSEVAL_ALT)]  # replace NA with lag
db1.MHSEVAL_ALT_SE[np.isnan(db1.MHSEVAL_ALT_SE)] = db1.MHSEVAL_ALT_LAG_SE[
    np.isnan(db1.MHSEVAL_ALT_SE)]  # replace NA with lag
db1['POPDENS_SE'] = se_ratio(db1.POPDENS, db.SE_T02A_002,
                             db.ACS12_5yr_B01003001s,
                             0)  # these are nearly all zero since county pops tend to have 0 MOE

# Unit test for equivalency
equal_test(db1.QUNOCCHU, db1.QUNOCCHU_ALT)
equal_test(db1.QED12LES, db1.QED12LES_ALT)
equal_test(db1.QEXTRCT, db1.QEXTRCT_ALT)
equal_test(db1.QSERV, db1.QSERV_ALT)
equal_test(db1.QNOAUTO, db1.QNOAUTO_ALT)

# Add in the sample sizes
db1['sample_person'] = db.ACS12_5yr_B00001001
db1['sample_hu'] = db.ACS12_5yr_B00002001

# Save data if main
if __name__ == "__main__":
    db1.to_csv(os.path.join(path, 'sovi_inputs.csv'))
