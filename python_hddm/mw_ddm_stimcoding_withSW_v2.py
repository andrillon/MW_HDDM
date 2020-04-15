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

mydata_face = mydata[mydata.Task==1]
mydata_digit = mydata[mydata.Task==2]

# I fit the baseline and state models to both the digit and faces tasks separately here. 
# I save the states, traces and traceplots. 
# For the two state models I added in some code for the posterior plots by condition. 
# PPC is a little more involved and is better to do when we've decided on a final set of models
#----------------------------------------#
#ColNames=['W_Fz','W_Cz','W_Pz','W_Oz']
ColNames=mydata.columns
#CentralElec = [11, 33, 22 ,26]
for nE in range(10, 63):
    print('***** WORKING ON ELEC ' + ColNames[nE] + ' *****\n\n\n')
#    # Both - Model 2 (SW)
#    model_W = hddm.HDDMStimCoding(mydata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', ColNames[nE]], 'a': ColNames[nE], 'z' : ColNames[nE], 't': ColNames[nE]}, p_outlier=.05)
#    model_W.find_starting_values()
#    model_W.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Both/model_SW/model_W.db'), db='pickle')
#    modelname = 'Both/model_SW/model_' + ColNames[nE]
#    #model_W.save(os.path.join(modelpath,modelname))
#    
#    # Extract stats
#    model_W_stats = model_W.gen_stats()
#    statsname = 'Both/model_SW/model_' + ColNames[nE] + '_stats.csv'
#    model_W_stats.to_csv(os.path.join(modelpath,statsname))
#    
#    # Extract traceplots
#    model_W.plot_posteriors(save=True, path=os.path.join(modelpath,'Both/model_SW/Plots/Traceplots'))
#    
#    # Extract full chains, transform back z and get v bias
#    model_W_traces = model_W.get_group_traces()
#    model_W_traces['z(0)'] = np.exp(model_W_traces['z_trans(0)'])/(1+np.exp(model_W_traces['z_trans(0)']))
#    model_W_traces['z(1)'] = np.exp(model_W_traces['z_trans(1)'])/(1+np.exp(model_W_traces['z_trans(1)']))
#    model_W_traces['v_bias(0)'] = abs(model_W_traces['v(0.1)'])-abs(model_W_traces['v(0.0)'])
#    model_W_traces['v_bias(1)'] = abs(model_W_traces['v(1.1)'])-abs(model_W_traces['v(1.0)'])
#    tracename = 'Both/model_SW/model_' + ColNames[nE] + '_traces.csv'
#    model_W_traces.to_csv(os.path.join(modelpath,tracename))
    
    # Face - Model 2 (SW)
    model_W = hddm.HDDMStimCoding(mydata_face, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', ColNames[nE]], 'a': ColNames[nE], 'z' : ColNames[nE], 't': ColNames[nE]}, p_outlier=.05)
    model_W.find_starting_values()
    model_W.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Faces/model_SW/model_W.db'), db='pickle')
    modelname = 'Faces/model_SW/model_' + ColNames[nE]
    #model_W.save(os.path.join(modelpath,modelname))    
    # Extract stats
    model_W_stats = model_W.gen_stats()
    statsname = 'Faces/model_SW/model_' + ColNames[nE] + '_stats.csv'
    model_W_stats.to_csv(os.path.join(modelpath,statsname))    
    # Extract traceplots
    model_W.plot_posteriors(save=True, path=os.path.join(modelpath,'Faces/model_SW/Plots/Traceplots'))
    
    # Extract full chains, transform back z and get v bias
    model_W_traces = model_W.get_group_traces()
    model_W_traces['z(0)'] = np.exp(model_W_traces['z_trans(0)'])/(1+np.exp(model_W_traces['z_trans(0)']))
    model_W_traces['z(1)'] = np.exp(model_W_traces['z_trans(1)'])/(1+np.exp(model_W_traces['z_trans(1)']))
    model_W_traces['v_bias(0)'] = abs(model_W_traces['v(0.1)'])-abs(model_W_traces['v(0.0)'])
    model_W_traces['v_bias(1)'] = abs(model_W_traces['v(1.1)'])-abs(model_W_traces['v(1.0)'])
    tracename = 'Faces/model_SW/model_' + ColNames[nE] + '_traces.csv'
    model_W_traces.to_csv(os.path.join(modelpath,tracename))
    
for nE in range(10, 63):
    print('***** WORKING ON ELEC ' + ColNames[nE] + ' *****\n\n\n')   
    # Digit - Model 2 (SW)
    model_W = hddm.HDDMStimCoding(mydata_digit, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', ColNames[nE]], 'a': ColNames[nE], 'z' : ColNames[nE], 't': ColNames[nE]}, p_outlier=.05)
    model_W.find_starting_values()
    model_W.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Digits/model_SW/model_W.db'), db='pickle')
    modelname = 'Digits/model_SW/model_' + ColNames[nE]
    #model_W.save(os.path.join(modelpath,modelname))
    
    # Extract stats
    model_W_stats = model_W.gen_stats()
    statsname = 'Digits/model_SW/model_' + ColNames[nE] + '_stats.csv'
    model_W_stats.to_csv(os.path.join(modelpath,statsname))    
    
    # Extract traceplots
    model_W.plot_posteriors(save=True, path=os.path.join(modelpath,'Digits/model_SW/Plots/Traceplots'))
    
    # Extract full chains, transform back z and get v bias
    model_W_traces = model_W.get_group_traces()
    model_W_traces['z(0)'] = np.exp(model_W_traces['z_trans(0)'])/(1+np.exp(model_W_traces['z_trans(0)']))
    model_W_traces['z(1)'] = np.exp(model_W_traces['z_trans(1)'])/(1+np.exp(model_W_traces['z_trans(1)']))
    model_W_traces['v_bias(0)'] = abs(model_W_traces['v(0.1)'])-abs(model_W_traces['v(0.0)'])
    model_W_traces['v_bias(1)'] = abs(model_W_traces['v(1.1)'])-abs(model_W_traces['v(1.0)'])
    tracename = 'Digits/model_SW/model_' + ColNames[nE] + '_traces.csv'
    model_W_traces.to_csv(os.path.join(modelpath,tracename))

