#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File for setting up plot labels and variables as dictionaries
Ensures order is consistent across plots
"""
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs ={'rcmf':False,'eta0':True,'beta':False,'etal':True,'Xs_0':False,'Fe0':True,'alpha_n':False,'r':False}
variables = ['etal','beta','alpha_n','rcmf','eta0','Xs_0','Fe0','r']
Myr = 365*24*3600*1e6 #number of s in Myr