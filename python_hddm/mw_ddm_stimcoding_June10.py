# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:07:24 2019

@author: Angus C. Burns

Code for generating a grid of variables based on slow-wave amplplitude and frequency 
and then fitting HDDMStimCoding models and saving output. 

You should be able to run the below code to preprocess the data, define the functions 
for fitting and model and then run the For loop at the bottom to run the grid search. 

The function will create a new folder for each of the models (currently 10 x 10 = 100) with 
a unique name based on the frequency F and amplitude A, e.g. model_splitF10A20.

Make sure to adjust the necessary paths:
    datapath
    datafile
    mypath 

Then you should be all set :-) !
    
The function currently saves the following:
    1. Traces for key parameters
    2. Statistics for the model parameters
    3. Model comparison indices DIC, AIC and BIC
    
You can also adjust the function in the for loop to run a longer chain if you like but 2k samples
seems sufficient from all the models I've fit so far. 
"""

## Import required packages
import hddm
import numpy as np         # for basic matrix operations
import pandas as pd
import os

# Load data from csv file into a NumPy structured array
datapath = 'D:\SlowWave\\'
datafile = 'HDDM_WIM_localsleep_pup_June5.txt'
data = mydata = pd.read_csv(os.path.join(datapath, datafile))

#Recode StimCat to cohere with upper and lower response boundaries
mydata.StimCat[mydata.StimCat == 0] = 2
mydata.StimCat[mydata.StimCat == 1] = 0
mydata.StimCat[mydata.StimCat == 2] = 1

# Fix column labels for HDDM
mydata = mydata.rename(index=str, columns ={"StimCat": "stim", "RT": "rt", "RCode": "response", "SubID": "subj_idx"})

#Pre-processing to remove RTs shorter than NDT (~.2s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .2]
mydata.rt[mydata.response == 0] = -1
mydata= mydata[mydata.rt < 1]
mydata.rt[mydata.response == 0] = 999

# Calculate total number of outliers 
total_outliers = len(data)-len(mydata)

# Create a series of new variables for splitting data based on freq and amp
freq = range(10,101,10)
amp = range(10,101,10)

for f in freq:
    for a in amp:
        mydata.loc[(mydata['pW_Fz'] >= a) & (mydata['pF_Fz'] >= f), 'splitF{0}A{1}'.format(f,a)] = 1
        mydata.loc[(mydata['pW_Fz'] < a) | (mydata['pF_Fz'] < f), 'splitF{0}A{1}'.format(f,a)] = 0  
        mydata.loc[(np.isnan(mydata['pW_Fz'])) | (np.isnan(mydata['pF_Fz'])), 'splitF{0}A{1}'.format(f,a)] = 0  

# Create a list for the model variables
splits = list(mydata.columns[28:128])

# Define functions to calculate aic and bic as additional model comparison metrics
def aic(self):
	k = len(self.get_stochastics())
	logp = sum([x.logp for x in self.get_observeds()['node']])	
	return 2 * k - 2 * logp

def bic(self):
    k = len(self.get_stochastics())
    n = len(self.data)
    logp = sum([x.logp for x in self.get_observeds()['node']])
    return -2 * logp + k * np.log(n)


    # ===================================================== #
    # Run_Model Function to create, fit & save HDDM Models
    # ===================================================== #

def run_model(mypath, model_name, split, n_samples, burn, thin):

    thispath = os.path.join(mypath, 'model_{}'.format(split))
    if not os.path.exists(thispath):
        os.mkdir(thispath)
        

    model_filename  = os.path.join(mypath, 'model_{}'.format(split), 'modelfit-{}.model'.format(split))
    print(model_filename)

    if model_name == 'stimcoding_base':

        # Model without pupil or slow-wave information for comparison
        m = hddm.HDDMStimCoding(mydata, stim_col='stim', split_param='z', p_outlier=0.05,
                include=('z'), depends_on={'v': ['stim']})
    elif model_name == 'stimcoding_z_SW': 
        
        m = hddm.HDDMStimCoding(mydata, stim_col='stim', split_param='z', p_outlier=0.05,
                include=('z'), depends_on={'v': ['stim', split], 'a': split, 't': split, 'z': split})
        
    print("finding starting values")
    m.find_starting_values() # this should help the sampling

    print("begin sampling")
    m.sample(n_samples, burn=burn, thin=thin, db='pickle', dbname=os.path.join(mypath, 'model_{}'.format(split), 'modelfit-{}.db'.format(split)))
    m.save(os.path.join(mypath, 'model_{}'.format(split), 'modelfit-{}.model'.format(split))) # save the model to disk

    # ============================================ #
    # save the output values
    # ============================================ #

    # save the model comparison indices
    df = dict()
    df['dic_original'] = [m.dic]
    df['aic'] = [aic(m)]
    df['bic'] = [bic(m)]
    df2 = pd.DataFrame(df)
    df2.to_csv(os.path.join(mypath, 'model_{}'.format(split), 'model_comparison_{}.csv'.format(split)))
    
    # save the stats
    results = m.gen_stats()
    results.to_csv(os.path.join(mypath, 'model_{}'.format(split), 'model_stats_{}.csv'.format(split)))    
    
    # save the traces for inspection
    params_of_interest = ['a(0.0)', 'a(1.0)', 'v(0.0.0.0)', 'v(1.0.0.0)', 'v(0.0.1.0)', 'v(1.0.1.0)', 't(0.0)', 't(1.0)', 'z(0.0)', 'z(1.0)']
    traces = []
    for p in range(len(params_of_interest)):
        traces.append(m.nodes_db.node[params_of_interest[p]].trace.gettrace())
    tracesarray = np.asarray(traces)
    tracesFrame= pd.DataFrame(data=tracesarray[0:,0:]) 
    tracesFrame.to_csv(os.path.join(mypath, 'model_{}'.format(split), 'traces_{}.csv'.format(split)))
    
    "---------------------------------------------------------------------------------------------"
    
    
# =============================================== #
# For loop to run the grid search and save output
# =============================================== #

mypath = 'D:\SlowWave\models\\'
model_name = 'stimcoding_z_SW'
    
modelCount = 0
for split in splits:
    modelCount = modelCount + 1
    print('We are up to model ', modelCount,'!')
    try:
        run_model(mypath, split,2000,1000,1)
    except:
        continue
    