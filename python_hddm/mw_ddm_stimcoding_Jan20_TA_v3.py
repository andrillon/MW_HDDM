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
modelpath = "/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models_TA/"
datafile = 'HDDM_WIM_localsleep_amp_pup_thrE90P2P_Dec21_v5.txt'
data = mydata = pd.read_csv(os.path.join(datapath, datafile))

# Fix column labels for HDDM
mydata = mydata.rename(index=str, columns ={"stimulus": "stim", "RT": "rt", "SubID": "subj_idx"})

#Pre-processing to remove RTs shorter than NDT (~.25s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .3]
mydata= mydata[mydata.DistProbe > -20]
#mydata.rt[mydata.response == 0] = -1
#mydata= mydata[mydata.rt < 1]
#mydata.rt[mydata.response == 0] = 999

# Calculate total number of outliers - 3143/82364 = 3.8%
total_outliers = len(data)-len(mydata)
print(total_outliers)

# Create separate datasets for the face (1) and digit (2) tasks
mydata_face = mydata[mydata.Task==1]
mydata_digit = mydata[mydata.Task==2]

# I fit the baseline and state models to both the digit and faces tasks separately here. 
# I save the states, traces and traceplots. 
# For the two state models I added in some code for the posterior plots by condition. 
# PPC is a little more involved and is better to do when we've decided on a final set of models

#----------------------------------------#

## Faces - Model 1 (baseline) 
#model_face_1 = hddm.HDDMStimCoding(mydata_face, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': 'stim'}, p_outlier=.05)
#model_face_1.find_starting_values()# Create model and start MCMC sampling
#model_face_1.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Faces/Model_1/model_face_1.db'), db='pickle')
#model_face_1.save(os.path.join(modelpath,'Faces/Model_1/model_face_1'))
#
## Extract stats
#model_face_1_stats = model_face_1.gen_stats()
#model_face_1_stats.to_csv(os.path.join(modelpath,'Faces/Model_1/model_face_1_stats.csv'))
#
## Extract traceplots
#model_face_1.plot_posteriors(save=True, path=os.path.join(modelpath,'Faces/Model_1/Plots/Traceplots'))
#
## Extract full chains, transform back z and get v bias
#model_face_1_traces = model_face_1.get_group_traces()
#model_face_1_traces['z'] = np.exp(model_face_1_traces['z_trans'])/(1+np.exp(model_face_1_traces['z_trans']))
#model_face_1_traces['v_bias'] = abs(model_face_1_traces['v(1)'])-abs(model_face_1_traces['v(0)'])
#model_face_1_traces.to_csv(os.path.join(modelpath,'Faces/Model_1/model_face_1_traces.csv'))
#

#----------------------------------------#

# Faces - Model 2 (state)
model_face_2 = hddm.HDDMStimCoding(mydata_face, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'State'], 'a': 'State', 'z' : 'State', 't': 'State'}, p_outlier=.05)
model_face_2.find_starting_values()
model_face_2.sample(100, burn=50, dbname=os.path.join(modelpath,'Faces/Model_2/model_face_2.db'), db='pickle')
#model_face_2.save(os.path.join(modelpath,'Faces/Model_2/model_face_2'))

# Extract stats
model_face_2_stats = model_face_2.gen_stats()
model_face_2_stats.to_csv(os.path.join(modelpath,'Faces/Model_2/model_face_2_stats.csv'))

# Extract traceplots
model_face_2.plot_posteriors(save=True, path=os.path.join(modelpath,'Faces/Model_2/Plots/Traceplots'))

# Extract full chains, transform back z and get v bias
model_face_2_traces = model_face_2.get_group_traces()
model_face_2_traces['z(1)'] = np.exp(model_face_2_traces['z_trans(1)'])/(1+np.exp(model_face_2_traces['z_trans(1)']))
model_face_2_traces['z(2)'] = np.exp(model_face_2_traces['z_trans(2)'])/(1+np.exp(model_face_2_traces['z_trans(2)']))
model_face_2_traces['z(3)'] = np.exp(model_face_2_traces['z_trans(3)'])/(1+np.exp(model_face_2_traces['z_trans(3)']))
model_face_2_traces['v_bias(1)'] = abs(model_face_2_traces['v(1.1)'])-abs(model_face_2_traces['v(1.0)'])
model_face_2_traces['v_bias(2)'] = abs(model_face_2_traces['v(2.1)'])-abs(model_face_2_traces['v(2.0)'])
model_face_2_traces['v_bias(3)'] = abs(model_face_2_traces['v(3.1)'])-abs(model_face_2_traces['v(3.0)'])
model_face_2_traces.to_csv(os.path.join(modelpath,'Faces/Model_2/model_faces_2_traces.csv'))

#----------------------------------------#

## Digits - Model 1 (baseline)  
#model_digit_1 = hddm.HDDMStimCoding(mydata_digit, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': 'stim'}, p_outlier=.05)
#model_digit_1.find_starting_values()# Create model and start MCMC sampling
#model_digit_1.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Digits/Model_1/model_digit_1_traces.db'), db='pickle')
##model_digit_1.save(os.path.join(modelpath,'Digits/Model_1/model_digit_1'))
#
## Extract stats
#model_digit_1_stats = model_digit_1.gen_stats()
#model_digit_1_stats.to_csv(os.path.join(modelpath,'Digits/Model_1/model_digit_1_stats.csv'))
#
## Extract traceplots
#model_digit_1.plot_posteriors(save=True, path=os.path.join(modelpath,'Digits/Model_1/Plots/Traceplots'))
#
## Extract full chains, transform back z and get v bias
#model_digit_1_traces = model_digit_1.get_group_traces()
#model_digit_1_traces['z'] = np.exp(model_digit_1_traces['z_trans'])/(1+np.exp(model_digit_1_traces['z_trans']))
#model_digit_1_traces['v_bias'] = abs(model_digit_1_traces['v(1)'])-abs(model_digit_1_traces['v(0)'])
#model_digit_1_traces.to_csv(os.path.join(modelpath,'Digits/Model_1/model_digit_1_traces.csv'))

#----------------------------------------#
# Digits - Model 2 (state) 
model_digit_2 = hddm.HDDMStimCoding(mydata_digit, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'State'], 'a': 'State', 'z' : 'State', 't': 'State'}, p_outlier=.05)
model_digit_2.sample(4000, burn=1000, dbname=os.path.join(modelpath,'Digits/Model_2/model_digit_2_traces.db'), db='pickle')
#model_digit_2.save(os.path.join(modelpath,'Digits/Model_2/model_digit_2'))

# Extract stats
model_digit_2_stats = model_digit_2.gen_stats()
model_digit_2_stats.to_csv(os.path.join(modelpath,'Digits/Model_2/model_digit_2_stats.csv'))

# Extract traceplots
model_digit_2.plot_posteriors(save=True, path=os.path.join(modelpath,'Digits/Model_2/Plots/Traceplots'))

# Extract full chains, transform z back and calculate v bias
model_digit_2_traces = model_digit_2.get_group_traces()
model_digit_2_traces['z(1)'] = np.exp(model_digit_2_traces['z_trans(1)'])/(1+np.exp(model_digit_2_traces['z_trans(1)']))
model_digit_2_traces['z(2)'] = np.exp(model_digit_2_traces['z_trans(2)'])/(1+np.exp(model_digit_2_traces['z_trans(2)']))
model_digit_2_traces['z(3)'] = np.exp(model_digit_2_traces['z_trans(3)'])/(1+np.exp(model_digit_2_traces['z_trans(3)']))
model_digit_2_traces['v_bias(1)'] = abs(model_digit_2_traces['v(1.1)'])-abs(model_digit_2_traces['v(1.0)'])
model_digit_2_traces['v_bias(2)'] = abs(model_digit_2_traces['v(2.1)'])-abs(model_digit_2_traces['v(2.0)'])
model_digit_2_traces['v_bias(3)'] = abs(model_digit_2_traces['v(3.1)'])-abs(model_digit_2_traces['v(3.0)'])
model_digit_2_traces.to_csv(os.path.join(modelpath,'Digits/Model_2/model_digit_2_traces.csv'))



