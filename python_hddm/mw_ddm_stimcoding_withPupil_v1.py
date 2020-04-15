# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:07:24 2019

@author: Angus C. Burns

Code for fitting the HDDM to 
"""

## Import required packages
import hddm
import numpy as np         # for basic matrix operations
import pandas as pd
import os
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt

# Load data from csv file into a NumPy structured array
datapath = "/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Data/"
modelpath = "/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/"
datafile = 'HDDM_WIM_localsleep_amp_pup_thrE90P2P_Dec21_v5.txt'
data = mydata = pd.read_csv(os.path.join(datapath, datafile))

# Fix column labels for HDDM
mydata = mydata.rename(index=str, columns ={"stimulus": "stim", "RT": "rt", "SubID": "subj_idx"})

#Pre-processing to remove RTs shorter than NDT (~.25s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .3]
#mydata= mydata[mydata.DistProbe > -19]
#mydata.rt[mydata.response == 0] = -1
#mydata= mydata[mydata.rt < 1]
#mydata.rt[mydata.response == 0] = 999

# Calculate total number of outliers - 3143/82364 = 3.8%
total_outliers = len(data)-len(mydata)
print(total_outliers)


# I fit the baseline and state models to both the digit and faces tasks separately here. 
# I save the states, traces and traceplots. 
# For the two state models I added in some code for the posterior plots by condition. 
# PPC is a little more involved and is better to do when we've decided on a final set of models
#----------------------------------------#

# Faces - Model 2 (state)
model_Pup = hddm.HDDMStimCoding(mydata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'pcPup'], 'a': 'pcPup', 'z' : 'pcPup', 't': 'pcPup'}, p_outlier=.05)
model_Pup.find_starting_values()
model_Pup.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Both/model_Pup/model_Pup.db'), db='pickle')
model_Pup.save(os.path.join(modelpath,'Both/model_Pup/model_Pup'))

# Extract stats
model_Pup_stats = model_Pup.gen_stats()
model_Pup_stats.to_csv(os.path.join(modelpath,'Both/model_Pup/model_Pup_stats.csv'))

# Extract traceplots
model_Pup.plot_posteriors(save=True, path=os.path.join(modelpath,'Both/model_Pup/Plots/Traceplots'))

# Extract full chains, transform back z and get v bias
model_Pup_traces = model_Pup.get_group_traces()
model_Pup_traces['z(1)'] = np.exp(model_Pup_traces['z_trans(1)'])/(1+np.exp(model_Pup_traces['z_trans(1)']))
model_Pup_traces['z(2)'] = np.exp(model_Pup_traces['z_trans(2)'])/(1+np.exp(model_Pup_traces['z_trans(2)']))
model_Pup_traces['z(3)'] = np.exp(model_Pup_traces['z_trans(3)'])/(1+np.exp(model_Pup_traces['z_trans(3)']))
model_Pup_traces['z(4)'] = np.exp(model_Pup_traces['z_trans(4)'])/(1+np.exp(model_Pup_traces['z_trans(4)']))
model_Pup_traces['z(5)'] = np.exp(model_Pup_traces['z_trans(5)'])/(1+np.exp(model_Pup_traces['z_trans(5)']))
model_Pup_traces['v_bias(1)'] = abs(model_Pup_traces['v(1.1)'])-abs(model_Pup_traces['v(1.0)'])
model_Pup_traces['v_bias(2)'] = abs(model_Pup_traces['v(2.1)'])-abs(model_Pup_traces['v(2.0)'])
model_Pup_traces['v_bias(3)'] = abs(model_Pup_traces['v(3.1)'])-abs(model_Pup_traces['v(3.0)'])
model_Pup_traces['v_bias(4)'] = abs(model_Pup_traces['v(4.1)'])-abs(model_Pup_traces['v(4.0)'])
model_Pup_traces['v_bias(5)'] = abs(model_Pup_traces['v(5.1)'])-abs(model_Pup_traces['v(5.0)'])
model_Pup_traces.to_csv(os.path.join(modelpath,'Both/model_Pup/model_Pup_traces.csv'))


