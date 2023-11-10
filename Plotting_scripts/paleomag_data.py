#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the paleomagnetic record
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Import data
paleo = pd.read_csv("../meteorite_paleomagnetism.csv",skiprows=[1])

#%% Process data
mclass = paleo['Classification'].unique()
yplot = np.arange(len(mclass)+1) #value to plot a class against

#plot all classes
for i, met in enumerate(mclass):
    mdata = paleo.loc[paleo['Classification']==met,:] #filter by class
    for j in range(len(mdata)): #plot each value in class - probably want this to be fill between 
    #also colour code by radio == True
        plt.scatter(mdata.loc[j,'rel_age_lower'],yplot[i])
        plt.scatter(mdata.loc[j,'rel_age_upper'],yplot[i])
    break

#Things to fix
# nan values
# what do we want on the plot
