#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for comparing latent heat and gravitational terms once have run solver
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#scale time to Myr
from parameters import Myr

run=17 # run with all terms in it

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
Tc = npzfile['Tc'] #core temp in K
t = npzfile['t'] #time in s
f = npzfile['f'] #time in s

t_plot = t/Myr

from ql_func import Qlt
from qg_func import Qgt

ql = Qlt(Tc,f)
qg = Qgt(Tc,f)

plt.figure()
plt.semilogy(t_plot,ql/qg)
plt.xlabel('t/Myr')
plt.ylabel('Ql/Qg')
plt.savefig('Plots/Latent_vs_grav.png')

