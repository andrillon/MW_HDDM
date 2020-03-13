# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:07:24 2019

@author: Angus C. Burns

Code for fitting the Go/No-Go HDDM to Digit and Face Task data. 
"""

## Import required packages
import hddm
import numpy as np         # for basic matrix operations
import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt

# Load data from csv file into a NumPy structured array
datapath = 'D:\Projects\MW_HDDM\Data\\'
modelpath = "D:\Projects\MW_HDDM\Models\\"
datafile = 'WanderIM_ProbeResults2_MW_new.txt'
data = mydata = pd.read_csv(os.path.join(datapath, datafile))

# Fix column labels for HDDM
mydata = mydata.rename(index=str, columns ={"stimulus": "stim", "RT": "rt", "SubID": "subj_idx"})

#Pre-processing to remove RTs shorter than NDT (~.25s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .3]
#mydata.rt[mydata.response == 0] = -1
#mydata= mydata[mydata.rt < 1]
#mydata.rt[mydata.response == 0] = 999

# Calculate total number of outliers - 3143/82364 = 3.8%
total_outliers = len(data)-len(mydata)
print(total_outliers)

# Create separate datasets for the face (1) and digit (2) tasks
pdata = mydata[np.logical_or(np.logical_and(mydata.DistProbe > -19,mydata.TrCat ==0),np.logical_and(mydata.DistProbe > -3,mydata.TrCat ==1))]

""" I titled the new sets of models omni short for omnibus as they all include both tasks in the same model"""

# Model 1 - All trials, Stim & Task only
# Model 2 - All trials, Stim, State & Task 
# Model 3 - Proximal trials (only 20 trials prior to probe), Stim & Task only
# Model 4 - Proximal trials (only 20 trials prior to probe), Stim, State & Task 

#--------------------------------------------------------------------------------------------------------#

# Omni - Model 1 (No State information, just Stim and Task) 
model_omni_1 = hddm.HDDMStimCoding(mydata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'Task'], 'a': 'Task', 'z' : 'Task', 't': 'Task'}, p_outlier=.05)
model_omni_1.find_starting_values()# Create model and start MCMC sampling
model_omni_1.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Omni\Model_1\model_omni_1.db'), db='pickle')
model_omni_1.save(os.path.join(modelpath,'Omni\Model_1\model_omni_1'))

# Extract stats
model_omni_1_stats = model_omni_1.gen_stats()
model_omni_1_stats.to_csv(os.path.join(modelpath,'Omni\Model_1\model_omni_1_stats.csv'))

# Extract traceplots
model_omni_1.plot_posteriors(save=True, path=os.path.join(modelpath,'Omni\Model_1\Plots\Traceplots'))

# Extract full chains, transform back z and get v bias
model_omni_1_traces = model_omni_1.get_group_traces()
model_omni_1_traces['z(1)'] = np.exp(model_omni_1_traces['z_trans(1)'])/(1+np.exp(model_omni_1_traces['z_trans(1)']))
model_omni_1_traces['z(2)'] = np.exp(model_omni_1_traces['z_trans(2)'])/(1+np.exp(model_omni_1_traces['z_trans(2)']))
model_omni_1_traces.to_csv(os.path.join(modelpath,'Omni\Model_1\model_omni_1_traces.csv'))

#----------------------------------------#

# Omni - Model 2 (State, Stim and Task) 
model_omni_2 = hddm.HDDMStimCoding(mydata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'State', 'Task'], 'a': ['State', 'Task'], 'z' : ['State', 'Task'], 't': ['State', 'Task']}, p_outlier=.05)
model_omni_2.find_starting_values()# Create model and start MCMC sampling
model_omni_2.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Omni\Model_2\model_omni_2.db'), db='pickle')
model_omni_2.save(os.path.join(modelpath,'Omni\Model_2\model_omni_2'))

# Extract stats
model_omni_2_stats = model_omni_2.gen_stats()
model_omni_2_stats.to_csv(os.path.join(modelpath,'Omni\Model_2\model_omni_2_stats.csv'))

# Extract traceplots
model_omni_2.plot_posteriors(save=True, path=os.path.join(modelpath,'Omni\Model_2\Plots\Traceplots'))

# Extract full chains, transform back z and get v bias
model_omni_2_traces = model_omni_2.get_group_traces()
model_omni_2_traces['z(1.1)'] = np.exp(model_omni_2_traces['z_trans(1.1)'])/(1+np.exp(model_omni_2_traces['z_trans(1.1)']))
model_omni_2_traces['z(2.1)'] = np.exp(model_omni_2_traces['z_trans(2.1)'])/(1+np.exp(model_omni_2_traces['z_trans(2.1)']))
model_omni_2_traces['z(3.1)'] = np.exp(model_omni_2_traces['z_trans(3.1)'])/(1+np.exp(model_omni_2_traces['z_trans(3.1)']))
model_omni_2_traces['z(1.2)'] = np.exp(model_omni_2_traces['z_trans(1.2)'])/(1+np.exp(model_omni_2_traces['z_trans(1.2)']))
model_omni_2_traces['z(2.2)'] = np.exp(model_omni_2_traces['z_trans(2.2)'])/(1+np.exp(model_omni_2_traces['z_trans(2.2)']))
model_omni_2_traces['z(3.2)'] = np.exp(model_omni_2_traces['z_trans(3.2)'])/(1+np.exp(model_omni_2_traces['z_trans(3.2)']))
model_omni_2_traces.to_csv(os.path.join(modelpath,'Omni\Model_2\model_omni_2_traces.csv'))

#----------------------------------------#

# Omni - Model 3 (state & task)
model_omni_3 = hddm.HDDMStimCoding(pdata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'Task'], 'a': 'Task', 'z' : 'Task', 't': 'Task'}, p_outlier=.05)
model_omni_3.find_starting_values()
model_omni_3.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Omni\Model_3\model_omni_3.db'), db='pickle')
model_omni_3.save(os.path.join(modelpath,'Omni\Model_3\model_omni_3'))
model_omni_3.print_stats()
model_omni_3.plot_posteriors_conditions()

# Extract stats
model_omni_3_stats = model_omni_3.gen_stats()
model_omni_3_stats.to_csv(os.path.join(modelpath,'Omni\Model_3\model_omni_3_stats.csv'))

# Extract traceplots
model_omni_3.plot_posteriors(save=True, path=os.path.join(modelpath,'Omni\Model_3\Plots\Traceplots'))

# Extract full chains, transform back z and get v bias
model_omni_3_traces = model_omni_3.get_group_traces()
model_omni_3_traces['z(1)'] = np.exp(model_omni_3_traces['z_trans(1)'])/(1+np.exp(model_omni_3_traces['z_trans(1)']))
model_omni_3_traces['z(2)'] = np.exp(model_omni_3_traces['z_trans(2)'])/(1+np.exp(model_omni_3_traces['z_trans(2)']))
model_omni_3_traces.to_csv(os.path.join(modelpath,'Omni\Model_3\model_omni_3_traces.csv'))

#----------------------------------------#

model_omni_4 = hddm.HDDMStimCoding(pdata, include={'z'}, stim_col= 'stim', split_param='z', depends_on={'v': ['stim', 'State', 'Task'], 'a': ['State', 'Task'], 'z' : ['State', 'Task'], 't': ['State', 'Task']}, p_outlier=.05)
model_omni_4.find_starting_values()
model_omni_4.sample(2000, burn=1000, dbname=os.path.join(modelpath,'Omni\Model_4\model_omni_4.db'), db='pickle')
model_omni_4.save(os.path.join(modelpath,'Omni\Model_4\model_omni_4'))
model_omni_4.print_stats()
model_omni_4.plot_posteriors_conditions()

# Extract stats
model_omni_4_stats = model_omni_4.gen_stats()
model_omni_4_stats.to_csv(os.path.join(modelpath,'Omni\Model_4\model_omni_4_stats.csv'))

# Extract traceplots
model_omni_4.plot_posteriors(save=True, path=os.path.join(modelpath,'Omni\Model_4\Plots\Traceplots'))

# Extract full chains, transform back z and get v bias
model_omni_4_traces = model_omni_4.get_group_traces()
model_omni_4_traces['z(1.1)'] = np.exp(model_omni_4_traces['z_trans(1.1)'])/(1+np.exp(model_omni_4_traces['z_trans(1.1)']))
model_omni_4_traces['z(2.1)'] = np.exp(model_omni_4_traces['z_trans(2.1)'])/(1+np.exp(model_omni_4_traces['z_trans(2.1)']))
model_omni_4_traces['z(3.1)'] = np.exp(model_omni_4_traces['z_trans(3.1)'])/(1+np.exp(model_omni_4_traces['z_trans(3.1)']))
model_omni_4_traces['z(1.2)'] = np.exp(model_omni_4_traces['z_trans(1.2)'])/(1+np.exp(model_omni_4_traces['z_trans(1.2)']))
model_omni_4_traces['z(2.2)'] = np.exp(model_omni_4_traces['z_trans(2.2)'])/(1+np.exp(model_omni_4_traces['z_trans(2.2)']))
model_omni_4_traces['z(3.2)'] = np.exp(model_omni_4_traces['z_trans(3.2)'])/(1+np.exp(model_omni_4_traces['z_trans(3.2)']))
model_omni_4_traces.to_csv(os.path.join(modelpath,'Omni\Model_4\model_omni_4_traces.csv'))

#----------------------------------------#


#---------------------------------------------------------------------------------------#

### Plots for State*Task Effects on Parameters - using Model 4

# Threshold
fig, axes = plt.subplots()
sns.kdeplot(model_omni_4_traces['a(1.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['a(2.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['a(3.1)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['a(1.2)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['a(2.2)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['a(3.2)'], vertical=False, shade=True,  ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
axes.set_xlim(xmin=1.2,xmax=2.9)
axes.set_title("State Effect on Decision Threshold", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Omni\Model_4\Plots\Threshold.jpeg'), dpi=600)

# Non-Decision Time
fig, axes = plt.subplots()
sns.kdeplot(model_omni_4_traces['t(1.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['t(2.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['t(3.1)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['t(1.2)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['t(2.2)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['t(3.2)'], vertical=False, shade=True,  ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
axes.set_xlim(xmin=0.14,xmax=.32)
axes.set_title("State Effect on NDT", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Omni\Model_4\Plots\\NDT.jpeg'), dpi=600)

# No-Go Drift
fig, axes = plt.subplots()
sns.kdeplot(model_omni_4_traces['v(1.1.0)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(2.1.0)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['v(3.1.0)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(1.2.0)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(2.2.0)'], vertical=False, shade=True, ax=axes)
sns.kdeplot(model_omni_4_traces['v(3.2.0)'], vertical=False, shade=True,  ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
axes.set_xlim(xmin=-2.5,xmax=0)
axes.set_title("State Effect on No-Go Drift", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Omni\Model_4\Plots\Drift_nogo.jpeg'), dpi=600)

# Go Drift
fig, axes = plt.subplots()
sns.kdeplot(model_omni_4_traces['v(1.1.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(2.1.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(3.1.1)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['v(1.2.1)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(2.2.1)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['v(3.2.1)'], vertical=False, shade=True, ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
axes.set_xlim(xmin=3.,xmax=6.)
axes.set_title("State Effect on Go Drift", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Omni\Model_4\Plots\Drift_go.jpeg'), dpi=600)

# Starting Point
fig, axes = plt.subplots()
sns.kdeplot(model_omni_4_traces['z(1.1)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['z(2.1)'], vertical=False, shade=False,  ax=axes)
sns.kdeplot(model_omni_4_traces['z(3.1)'], vertical=False, shade=False, ax=axes)
sns.kdeplot(model_omni_4_traces['z(1.2)'], vertical=False, shade=True, ax=axes)
sns.kdeplot(model_omni_4_traces['z(2.2)'], vertical=False, shade=True,  ax=axes)
sns.kdeplot(model_omni_4_traces['z(3.2)'], vertical=False, shade=True,  ax=axes)
axes.set_xlabel('Parameter Estimate a.u.', fontsize='large', fontweight='bold')
axes.set_xlim(xmin=0.16,xmax=.42)
axes.set_title("State Effect on Starting Point", fontsize='x-large', fontweight='bold')
axes.set_ylabel('Posterior Probability', fontsize='large', fontweight='bold')
axes.tick_params(axis='both',direction='out',length=6,width=2, labelsize='large')
sns.despine()
fig.savefig(os.path.join(modelpath, 'Omni\Model_4\Plots\StartingPoint.jpeg'), dpi=600)

# Example code for calculating overlap in posterior density
a_ON, a_MW, a_MB = model_omni_2.nodes_db.node[['a(1)', 'a(2)' , 'a(3)']]
print("P_v(ON > MB) = ", (a_ON.trace() > a_MW.trace()).mean()) # P = 0.98