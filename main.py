#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for executing multiple runs
"""
import runpy
import importlib
import parameters
nruns = 6 #how many runs are there

i=0
while i < nruns:
    runpy.run_path('solver.py')
    #reload parameters file
    importlib.reload(parameters)
    i +=1
    print(f'Run {i} complete')
    