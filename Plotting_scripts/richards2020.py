#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contour plot based on Fig 7 in Richards et. al. (2020)
"""
#currently a bit sparse as will need more data
import matplotlib.pyplot as plt
import pandas as pd

from filter_data import filter_two
from scatter_function import make_sub_scatter

folder = 'Fullrun2'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km

varlist =['r','Xs_0','rcmf','frht','etal','eta0','Fe0']
varlabel =['radius /km','$X_{S,0}$ /wt%','rcmf','frht','liquid viscosity /Pas','reference viscosity /Pas','$^{60}Fe/^{56}Fe$']
logvar = [False,False,False,False,True,True,False]

surf_var = 'tsolid' #what is the surface plot
surf_label = 'Solidification time /Myr'
surfmin = min(data[surf_var])
surfmax = max(data[surf_var])
save = False #do you want to save the figures
#fixed values for when quantities vary
eta0=1e21
r = 400
rcmf = 0.3
frht=0.005
Xs_0 = 28.5
Fe0 = 0
alpha_n = 25
etal = 10

i=0
j=0
var2 = 'r'
varlab2 = 'radius /km'
fig, ax = plt.subplots(nrows=len(varlist)-1,ncols=len(varlist)-1,tight_layout=True,figsize=[10,10],sharex='col',sharey='row')
for var2, lab2, log2 in zip(varlist[:-1],varlabel[:-1],logvar[:-1]):
    for var1, lab1, log1 in zip(varlist[j+1:],varlabel[j+1:],logvar[j+1:]):
            data_fil = filter_two(data,var1,var2,r,Xs_0,Fe0,rcmf,frht,etal,eta0,alpha_n)
            datap = data_fil.loc[:,[var1,var2,surf_var]] #make subset of data
            datap = datap.pivot(index=var1,columns=var2,values=surf_var)
            #make plot
            p1=ax[i,j].pcolormesh(datap.columns,datap.index,datap,vmin=surfmin,vmax=surfmax)
            cs=ax[i,j].contour(datap.columns,datap.index,datap,4,linestyles='--',colors='white',linewidths=0.5)
            plt.clabel(cs,fmt='%2.0f')
            if log2 == True:
                ax[i,j].set_xscale('log')
            if var1 == 'etal': #fix liquid viscosity axis scale bug
                ax[i,j].set_ylim([10,1000])
            if log1 == True:
                ax[i,j].set_yscale('log')
            if j == 0:
                ax[i,j].set_ylabel(lab1)
            if i==len(varlist)-2:
                ax[i,j].set_xlabel(lab2)
            i=i+1
    j=j+1
    i = j
fig.colorbar(p1,ax=ax[:,len(varlist)-2],label=surf_label)
fig.suptitle(f'Effect on {surf_label} \n Constant values: r={r} km, $X_{{S,0}}$={Xs_0} wt%, $^{{60}}Fe/^{{56}}Fe$={Fe0} \n $\\phi_{{rcmf}}$={rcmf},frht={frht} $K^{{-1}}$, $\\eta_l$={etal} Pas, $\\eta_0$={eta0} Pas')
#remove unused axes
for i in range(len(varlist)-1):
    for j in range(i+1,len(varlist)-1):
        ax[i,j].remove()
plt.savefig(f'../Plots/{folder}/multivar_surf_{surf_var}.png',dpi=450)


######################### Another option using contourf ########################
# i=0
# j=0
# var2 = 'r'
# varlab2 = 'radius /km'
# fig, ax = plt.subplots(nrows=len(varlist)-1,ncols=len(varlist)-1,tight_layout=True,figsize=[10,10],sharex='col',sharey='row')
# for var2, lab2, log2 in zip(varlist[:-1],varlabel[:-1],logvar[:-1]):
#     for var1, lab1, log1 in zip(varlist[j+1:],varlabel[j+1:],logvar[j+1:]):
#             data_fil = filter_two(data,var1,var2,r,Xs_0,Fe0,rcmf,frht,etal,eta0,alpha_n)
#             datap = data_fil.loc[:,[var1,var2,surf_var]] #make subset of data
#             datap = datap.pivot(index=var1,columns=var2,values=surf_var)
#             #make plot
#             if (i==0) & (j==0): #use first set for colourbar
#                 cs1=ax[i,j].contourf(datap.columns,datap.index,datap,4,vmin=surfmin,vmax=surfmax)
#                 plt.clabel(cs1,fmt='%2.0f',colors='white')
#             else:
#                 cs=ax[i,j].contourf(datap.columns,datap.index,datap,4,vmin=surfmin,vmax=surfmax)
#                 plt.clabel(cs,fmt='%2.0f',colors='white')
#             if log2 == True:
#                 ax[i,j].set_xscale('log')
#             if var1 == 'etal': #fix liquid viscosity axis scale bug
#                 ax[i,j].set_ylim([10,1000])
#             if log1 == True:
#                 ax[i,j].set_yscale('log')
#             if j == 0:
#                 ax[i,j].set_ylabel(lab1)
#             if i==len(varlist)-2:
#                 ax[i,j].set_xlabel(lab2)
#             i=i+1
#     j=j+1
#     i = j
# fig.colorbar(cs1,ax=ax[:,len(varlist)-2],label=surf_label)
# fig.suptitle(f'Effect on {surf_label} \n Constant values: r={r} km, $X_{{S,0}}$={Xs_0} wt%, $^{{60}}Fe/^{{56}}Fe$={Fe0} \n $\\phi_{{rcmf}}$={rcmf},frht={frht} $K^{{-1}}$, $\\eta_l$={etal} Pas, $\\eta_0$={eta0} Pas')
# #remove unused axes
# for i in range(len(varlist)-1):
#     for j in range(i+1,len(varlist)-1):
#         ax[i,j].remove()
        
################# Same grid but for a scatter plot ############################
i=0
j=0
var2 = 'r'
varlab2 = 'radius /km'
fig, ax = plt.subplots(nrows=len(varlist)-1,ncols=len(varlist)-1,tight_layout=True,figsize=[10,10],sharex='col',sharey='row')
for var2, lab2, log2 in zip(varlist[:-1],varlabel[:-1],logvar[:-1]):
    for var1, lab1, log1 in zip(varlist[j+1:],varlabel[j+1:],logvar[j+1:]):
            data_fil = filter_two(data,var1,var2,r,Xs_0,Fe0,rcmf,frht,etal,eta0,alpha_n)
            p1 = ax[i,j].scatter(data_fil[var2], data_fil[var1],c=data_fil[surf_var],vmin=surfmin,vmax=surfmax)
            if log2 == True:
                ax[i,j].set_xscale('log')
            if var1 == 'etal': #fix liquid viscosity axis scale bug
                ax[i,j].set_ylim([10,1000])
            if log1 == True:
                ax[i,j].set_yscale('log')
            if j == 0:
                ax[i,j].set_ylabel(lab1)
            if i==len(varlist)-2:
                ax[i,j].set_xlabel(lab2)
            i=i+1
    j=j+1
    i = j
fig.colorbar(p1,ax=ax[:,len(varlist)-2],label=surf_label)
fig.suptitle(f'Effect on {surf_label} \n Constant values: r={r} km, $X_{{S,0}}$={Xs_0} wt%, $^{{60}}Fe/^{{56}}Fe$={Fe0} \n $\\phi_{{rcmf}}$={rcmf},frht={frht} $K^{{-1}}$, $\\eta_l$={etal} Pas, $\\eta_0$={eta0} Pas')
#remove unused axes
for i in range(len(varlist)-1):
    for j in range(i+1,len(varlist)-1):
        ax[i,j].remove()

plt.savefig(f'../Plots/{folder}/multivar_scatter_{surf_var}.png',dpi=450)