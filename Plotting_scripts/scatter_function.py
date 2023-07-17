#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for creating a scatter plot with a 3rd dimension specified by colour and a fourth by size
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def make_scatter(x,y,xlabel,ylabel,colour='black',size=10,ss=10,log=[False, False, False], colourlabel=None,sizelabel=None,save=[False,'path']):
    """
    Draw a scatter plot

    Parameters
    ----------
    x : np.array
        x axis variable
    y : np.array
        y axis variable
    xlabel : string
        label for x axis
    ylabel : string
        label for y axis
    colour : np.array
        variable to represented by colour of points. The default is 'black'.
    size : float, optional
        varaible represented by size of points. The default is 5
    ss : float
        size scale, amount to multiply size scale by to get reasonably sized points
    log : list of bool, optional
        should x, y and colour map be log scaled [x, y, colormap]
    colourlabel : string optional
        label for colorbar The default is None.
    sizelabel : string, optional
        label for point size scale. The default is None.
    save : [bool, str]
        should the figure be saved and filename and path to save figure

    Returns
    -------
    None.

    """
    
    plt.figure()
    if log[2] == True:
        plt.scatter(x,y,s=ss*size,c=colour,norm=mcolors.LogNorm())
    else:
        plt.scatter(x,y,s=ss*size,c=colour)
    if log[0] == True:
        plt.xscale('log')
    if log[1] == True:
        plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colourlabel != None:
        plt.colorbar(label=colourlabel) 
    if sizelabel != None:
        for val in [min(size), np.mean(size), max(size)]:
            plt.scatter([], [], c='k', alpha=0.3, s=ss*val,label=f'{val:.2f}')
            plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title=sizelabel)
    if save[0] == True:
        plt.savefig(save[1],dpi=450)
    return None

def make_sub_scatter(x,y,xlabel,ylabel,row,column,n,colour='lightseagreen',size=10,ss=10,log=[False, False, False], colourlabel=None,sizelabel=None):
    """
    Draw a scatter plot as a subplot

    Parameters
    ----------
    x : np.array
        x axis variable
    y : np.array
        y axis variable
    xlabel : string
        label for x axis
    ylabel : string
        label for y axis
    row : int
        number of rows in subplots
    col : int 
        number of columns in subplots
    n : int
        position in subplots
    colour : np.array
        variable to represented by colour of points. The default is 'black'.
    size : float, optional
        variable represented by size of points. The default is 5
    ss : float
        size scale, amount to multiply size scale by to get reasonably sized points
    log : list of bool, optional
        should x, y and colour map be log scaled [x, y, colormap]
    colourlabel : string optional
        label for colorbar The default is None.
    sizelabel : string, optional
        label for point size scale. The default is None.

    Returns
    -------
    None.

    """
    
    plt.subplot(row,column,n)
    if log[2] == True:
        plt.scatter(x,y,s=ss*size,c=colour,norm=mcolors.LogNorm())
    else:
        plt.scatter(x,y,s=ss*size,c=colour)
    if log[0] == True:
        plt.xscale('log')
    if log[1] == True:
        plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colourlabel != None:
        plt.colorbar(label=colourlabel)
    if sizelabel != None:
        for val in [min(size), np.mean(size), max(size)]:
            plt.scatter([], [], c='k', alpha=0.3, s=ss*val,label=f'{val:.0f}')
            plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title=sizelabel)

    return None