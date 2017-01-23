import numpy as np
import os
import pandas as pd
from scipy.stats.mstats import zscore as ZSCORE
from spss_pca import SPSS_PCA

def getInputNames():

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
                   
    return input_names

    
def computeSoVI(data,input_names,outPath,region,attr,compute_sovi=True,output_weights=False,simulation=False,sims=None):
    
    """
    data: a pandas dataframe containing ACS vars written to counties
            input should have a county ID column named 'Geo_FIPS'
    
    input names: a nested list of input attributes, directionalities, and sample types
    
    outPath: output directory    
    
    region: study region - a string, for naming ouputs
    
    compute_sovi: if True (default) calculates SoVI values & writes output
    
    output_weights: if True writes the sum of weights (contribution to SoVI) for  
                    the given geography and variables. 
                    
    attr: attribute type used in analysis (i.e. all attributes vs. drop-one)
            - a string, for naming ouputs
    
    simulation: if True simulates SoVI over each variable [sims] times 
    
    sims: number of simulations
    """

    # some data frames to organize the input data
    inputs = pd.DataFrame(index=data.Geo_FIPS)
    zinputs = pd.DataFrame(index=data.Geo_FIPS)
    errors = pd.DataFrame(index=data.Geo_FIPS)
    
    # input data prep
    # --swap signs of the attributes expected to have a "negative" affect on vulnerability
    # --take z-score of data
    for name, sign, sample in input_names:
        if sign == 'neg':
            inputs[name] = -data[name].values
        elif sign == 'pos':
            inputs[name] = data[name].values
        else:
            raise Exception, "problem"
        errors[name] = data[name+'_SE'].values
        zinputs[name] = ZSCORE(inputs[name].values)
    zinputs_static = zinputs.copy(deep=True)
    
    if compute_sovi:    
    
        # compute SoVI
        pca = SPSS_PCA(zinputs, reduce=True, varimax=True)
        sovi_actual = pca.scores_rot.sum(1)
        sovi_actual = pd.DataFrame(sovi_actual, index=data.Geo_FIPS, columns=['sovi'])
        
        # write SoVI to hard drive
        sovi_actual.to_csv(os.path.join(os.getcwd(),outPath,'SoVI_'+region+'_'+attr+'.csv'))
        
    
    if output_weights:
        
        # share of each attribute that goes into SoVI
        # NOTE: since attributes are on different scales these shares cannot be
        #       compared across attributes
        attrib_contribution = pca.weights_rot.sum(1)
        sovi_alt_computation = (zinputs * attrib_contribution).sum(1) # this is just a check
        sovi_alt_computation = pd.DataFrame(sovi_alt_computation, columns=['sovi'])
        if not np.allclose(sovi_actual, sovi_alt_computation):
            raise Exception, "mismatch"
            
        ##Write attribute contribution output 
        
        #Get attribute names
        attr_names=[i[0] for i in input_names]
        
        #Generate dictionary
        c={'attribute':attr_names,'rot_sum':attrib_contribution}
        #Generate pandas dataframe
        contrib=pd.DataFrame(c)
        #Write to csv 
        outCSV=region+'_attribute+contrib_'+attr+'.csv'
#        contrib.to_csv(outPath+outCSV)
        contrib.to_csv(os.path.join(os.getcwd(),outPath,outCSV))    
        print "Wrote",outCSV+".\n"       
        
        ##Write factor loadings to CSV 
        
        loadings=pd.DataFrame(pca.weights_rot,index=attr_names)
        loadings.to_csv(os.path.join(os.getcwd(),outPath,outCSV+'_loadings.csv'))
    
    if simulation:
        
        ###########################################
        #### Compute simulated values of SoVI  ####
        ###########################################
        
        np.random.seed(87654)
        # take random draws on one attribute at a time
        for sim in range(sims):
            #print "Running",str(sims),"simulations of",name,"...\n"
            # some data frames to organize inputs and outputs
            colnames = [aName for aName, sign, sample in input_names]
            sim_attribs = pd.DataFrame(index=data.Geo_FIPS, columns=colnames)
            sim_sovis = pd.DataFrame(index=data.Geo_FIPS, columns=colnames)
            # get all simulated values for each county
            for name, sign, sample in input_names:
                for county in zinputs.index:
                    sample_size = data.ix[county]['sample_'+sample]
                    est = inputs.ix[county][name]
                    #From Brault (2014) that crosses from replicate weights SE to SD
                    err = np.sqrt((errors.ix[county][name]**2) * ((1-0.5)**2))  
                    if err > 0:   # need to catch the case of attributes with no uncertainty
                    # draw simulated value
                        sim_attribs.ix[county][name] = np.random.normal(est, err, sims)
                    else:
                    # return the value since no uncertainty
                        sim_attribs.ix[county][name] = est
            # prep to get data ready for PCA
                sim_attribs = sim_attribs.values.astype(float)
                sim_attribs = ZSCORE(sim_attribs)
                #zinputs = pd.DataFrame(sim_attribs, index=data.Geo_FIPS, columns=columns)
            # compute SoVI for each simulation
            pca = SPSS_PCA(sim_attribs, reduce=True, varimax=True)
            sovi = np.dot(zinputs_static, pca.weights_rot).sum(1)
            sim_sovis['sim'+str(sim)] = sovi
            # write results to hard drive
            #sim_sovis.to_csv(path+'results_counties/'+name+'_'+str(sims)+'.csv')
#            outSim=outPath+region+'_'+name+'_'+str(sims)+'.csv'
            outSim=os.path.join(os.getcwd(),outPath,region+'_'+name+'_'+str(sims)+'.csv')    
            sim_sovis.to_csv(outSim)
            print "Wrote",outSim+".\n"
            # cleanup for next attribute
            zinputs = zinputs_static.copy(deep=True)