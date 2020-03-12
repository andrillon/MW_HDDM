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
datafile = 'WanderIM_ProbeResults2_MW_new.txt'
data = mydata = pd.read_csv(os.path.join(datapath, datafile))

# Fix column labels for HDDM
mydata = mydata.rename(index=str, columns ={"stimulus": "stim", "RT": "rt", "SubID": "subj_idx"})

#Pre-processing to remove RTs shorter than NDT (~.25s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .3]
mydata= mydata[np.logical_or(np.logical_and(mydata.DistProbe > -19,mydata.TrCat ==0),np.logical_and(mydata.DistProbe > -3,mydata.TrCat ==1))]
#mydata.rt[mydata.response == 0] = -1
#mydata= mydata[mydata.rt < 1]
#mydata.rt[mydata.response == 0] = 999

# Calculate total number of outliers - 3143/82364 = 3.8%
total_outliers = len(data)-len(mydata)
print(total_outliers)

# Create separate datasets for the face (1) and both (2) tasks
mydata_both = mydata

# I fit the baseline and state models to both the both and faces tasks separately here. 
# I save the states, traces and traceplots. 
# For the two state models I added in some code for the posterior plots by condition. 
# PPC is a little more involved and is better to do when we've decided on a final set of models

#----------------------------------------#

# Digits - Model 1 (baseline)  
model_both_1 = hddm.HDDMStimCoding(mydata_both, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['Task', 'stim']}, p_outlier=.05)
model_both_1.find_starting_values()# Create model and start MCMC sampling
model_both_1.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Both/Model_1/model_both_1_traces.db'), db='pickle')
model_both_1.save(os.path.join(modelpath,'Both/Model_1/model_both_1'))

# Extract stats
model_both_1_stats = model_both_1.gen_stats()
model_both_1_stats.to_csv(os.path.join(modelpath,'Both/Model_1/model_both_1_stats.csv'))

# Extract traceplots
model_both_1.plot_posteriors(save=True, path=os.path.join(modelpath,'Both/Model_1/Plots/Traceplots'))

# Extract full chains, transform back z and get v bias
model_both_1_traces = model_both_1.get_group_traces()
model_both_1_traces['z'] = np.exp(model_both_1_traces['z_trans'])/(1+np.exp(model_both_1_traces['z_trans']))
model_both_1_traces['v_bias'] = abs(model_both_1_traces['v(1)'])-abs(model_both_1_traces['v(0)'])
model_both_1_traces.to_csv(os.path.join(modelpath,'Both/Model_1/model_both_1_traces.csv'))

#----------------------------------------#
# boths - Model 2 (state) 
model_both_2 = hddm.HDDMStimCoding(mydata_both, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['Task', 'stim', 'State'], 'a': ['Task', 'State'], 'z' : ['Task', 'State'], 't': ['Task', 'State']}, p_outlier=.05)
model_both_2.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Both/Model_2/model_both_2_traces.db'), db='pickle')
model_both_2.save(os.path.join(modelpath,'Both/Model_2/model_both_2'))

# Extract stats
model_both_2_stats = model_both_2.gen_stats()
model_both_2_stats.to_csv(os.path.join(modelpath,'Both/Model_2/model_both_2_stats.csv'))

# Extract traceplots
model_both_2.plot_posteriors(save=True, path=os.path.join(modelpath,'Both/Model_2/Plots/Traceplots'))

# Extract full chains, transform z back and calculate v bias
model_both_2_traces = model_both_2.get_group_traces()
model_both_2_traces['z(1.1)'] = np.exp(model_both_2_traces['z_trans(1.1)'])/(1+np.exp(model_both_2_traces['z_trans(1.1)']))
model_both_2_traces['z(2.1)'] = np.exp(model_both_2_traces['z_trans(2.1)'])/(1+np.exp(model_both_2_traces['z_trans(2.1)']))
model_both_2_traces['z(3.1)'] = np.exp(model_both_2_traces['z_trans(3.1)'])/(1+np.exp(model_both_2_traces['z_trans(3.1)']))
model_both_2_traces['z(1.2)'] = np.exp(model_both_2_traces['z_trans(1.2)'])/(1+np.exp(model_both_2_traces['z_trans(1.2)']))
model_both_2_traces['z(2.2)'] = np.exp(model_both_2_traces['z_trans(2.2)'])/(1+np.exp(model_both_2_traces['z_trans(2.2)']))
model_both_2_traces['z(3.2)'] = np.exp(model_both_2_traces['z_trans(3.2)'])/(1+np.exp(model_both_2_traces['z_trans(3.2)']))
model_both_2_traces['v_bias(1.1)'] = abs(model_both_2_traces['v(1.1.1)'])-abs(model_both_2_traces['v(1.1.0)'])
model_both_2_traces['v_bias(2.1)'] = abs(model_both_2_traces['v(2.1.1)'])-abs(model_both_2_traces['v(2.1.0)'])
model_both_2_traces['v_bias(3.1)'] = abs(model_both_2_traces['v(3.1.1)'])-abs(model_both_2_traces['v(3.1.0)'])
model_both_2_traces['v_bias(1.2)'] = abs(model_both_2_traces['v(1.2.1)'])-abs(model_both_2_traces['v(1.2.0)'])
model_both_2_traces['v_bias(2.2)'] = abs(model_both_2_traces['v(2.2.1)'])-abs(model_both_2_traces['v(2.2.0)'])
model_both_2_traces['v_bias(3.2)'] = abs(model_both_2_traces['v(3.2.1)'])-abs(model_both_2_traces['v(3.2.0)'])
model_both_2_traces.to_csv(os.path.join(modelpath,'Both/Model_2/model_both_2_traces.csv'))

#---------------------------------------------------------------------------------------#

### Plots for State Effects on Parameters for boths_2 & Faces_2 Models

###--- boths_2 Plots ---###

#Non-Decision Time
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['t(1.1)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['t(2.1)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['t(3.1)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['t(1.2)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['t(2.2)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['t(3.2)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=.22,xmax=.32)
axes.set_title("State Effect on NDT", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots//NDT.jpeg'), dpi=300)

# Threshold
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['a(1.1)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['a(2.1)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['a(3.1)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['a(1.2)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['a(2.2)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['a(3.2)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=1.2,xmax=2.2)
axes.set_title("State Effect on Decision Threshold", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots/Threshold.jpeg'), dpi=300)

# Starting Point
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['z(1.1)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['z(2.1)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['z(3.1)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['z(1.2)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['z(2.2)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['z(3.2)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=0.2,xmax=.45)
axes.set_title("State Effect on Starting Point", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots/StartingPoint.jpeg'), dpi=300)

# Drift Bias
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['v_bias(1.1)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v_bias(2.1)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v_bias(3.1)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['v_bias(1.2)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v_bias(2.2)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v_bias(3.2)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=2.5,xmax=4.8)
axes.set_title("State Effect on Drift Bias", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots/DriftBias.jpeg'), dpi=300)

# Go Drift
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['v(1.1.1)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v(2.1.1)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v(3.1.1)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['v(1.2.1)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v(2.2.1)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v(3.2.1)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=3.5,xmax=6.5)
axes.set_title("State Effect on Go Drift", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots/Drift_go.jpeg'), dpi=300)

# No-Go Drift
fig, axes = plt.subplots()
sns.kdeplot(model_both_2_traces['v(1.1.0)'], vertical=False, shade=False, color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v(2.1.0)'], vertical=False, shade=False, color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v(3.1.0)'], vertical=False, shade=False, color='b', ax=axes)
sns.kdeplot(model_both_2_traces['v(1.2.0)'], vertical=False, shade=False, linestyle="--", color='r', ax=axes)
sns.kdeplot(model_both_2_traces['v(2.2.0)'], vertical=False, shade=False, linestyle="--", color='g', ax=axes)
sns.kdeplot(model_both_2_traces['v(3.2.0)'], vertical=False, shade=False, linestyle="--", color='b', ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
#axes.set_xlim(xmin=-2.5,xmax=0)
axes.set_title("State Effect on No-Go Drift", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Both/Model_2/Plots/Drift_nogo.jpeg'), dpi=300)

