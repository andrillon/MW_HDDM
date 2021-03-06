# -*- coding: utf-8 -*-
"""
Created on Sun May 12 15:10:18 2019

@author: Angus C. Burns

Analysis of Mind Wandering Data for Slow-Wave and Pupil effects
"""

## Import required packages
import hddm
from patsy import dmatrix  # for generation of (regression) design matrices
import numpy as np         # for basic matrix operations
import pandas as pd
from pandas import Series  # to manipulate data-frames generated by hddm
import matplotlib.pyplot as plt
import os, sys, pickle, time
import seaborn as sns
import pylab as pl
from kabuki.analyze import gelman_rubin
from scipy import stats

# Load data from csv file into a NumPy structured array
data = mydata = hddm.load_csv('GNG_MW_DDM.csv')

#Pre-processing to remove RTs shorter than NDT (~.2s) and longer than stimulus duration (1s)
mydata.rt[mydata.response == 0] = 999
mydata= mydata[mydata.rt > .2]
mydata.rt[mydata.response == 0] = -1
mydata= mydata[mydata.rt < 1]

#Create new datasets for regression models 
Fz_mydata = mydata[mydata.W_Fz > 20]
Cz_mydata = mydata[mydata.W_Cz > 20]
Pz_mydata = mydata[mydata.W_Pz > 20]
Oz_mydata = mydata[mydata.W_Oz > 20]

#Z-score regressors
Fz_mydata['z_W_Fz'] = stats.zscore(Fz_mydata.W_Fz)
Cz_mydata['z_W_Cz'] = stats.zscore(Cz_mydata.W_Cz)
Pz_mydata['z_W_Pz'] = stats.zscore(Pz_mydata.W_Pz)
Oz_mydata['z_W_Oz'] = stats.zscore(Oz_mydata.W_Oz)


## Define link functions - Note: so far only v_link_func works

def z_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.stim.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * stim)))

def v_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stim.ix[x.index]})))
    return x * stim

#  Definte regression equations
v_reg = {'model': 'v ~ 1 + z_W_Pz:C(stim)', 'link_func': v_link_func}
a_reg = {'model': 'a ~ z_W_Pz', 'link_func': lambda a: a}
t_reg = {'model': 't ~ z_W_Pz', 'link_func': lambda t: t}
reg_descr = [v_reg, a_reg, t_reg]

# Fit model and print stats/posteriors
m_reg = hddm.HDDMRegressor(Pz_mydata, reg_descr, include='z')
m_reg.find_starting_values()
m_reg.sample(2000, burn=700)
m_reg.print_stats()
m_reg.plot_posteriors()
m_reg.plot_posteriors_conditions()

# Plot Posterior Nodes for Regressors, e.g. these plot for Fz nodes
SW_Go_v = m_reg.nodes_db.node["v_z_W_Fz:C(stim)[1]"]
hddm.analyze.plot_posterior_nodes([SW_Go_v], bins=20)
plt.xlabel('Slow-wave coefficient: Go Drift Rate ')
print("P(Go_v_theta < 0) = ", (SW_Go_v.trace() < 0).mean())

SW_NoGo_v = m_reg.nodes_db.node["v_z_W_Fz:C(stim)[0]"]
hddm.analyze.plot_posterior_nodes([SW_NoGo_v], bins=20)
plt.xlabel('Slow-Wave coeffecient: No-Go Drift Rate ')
print("P(NoGo_v_theta < 0) = ", (SW_NoGo_v.trace() < 0).mean())

SW_Boundary = m_reg.nodes_db.node["a_z_W_Fz"]
hddm.analyze.plot_posterior_nodes([SW_Boundary], bins=20)
plt.xlabel('Slow-Wave coeffecient: Boundary')
print("P(a_theta < 0) = ", (SW_Boundary.trace() < 0).mean())

SW_NDT = m_reg.nodes_db.node["t_z_W_Fz"]
hddm.analyze.plot_posterior_nodes([SW_NDT], bins=20)
plt.xlabel('Slow-Wave coeffecient: Non-Decision Time ')
print("P(t_theta < 0) = ", (SW_NDT.trace() < 0).mean())

