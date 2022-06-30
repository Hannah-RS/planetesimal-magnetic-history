#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for time evolution of CMB temperature and inner core radius based on Nimmo 2009
Uses the assumption Ephi = Eg not the full equations

Employs fourth order Runge-Kutta to solve Tc as a function of time

Creates 3 plots:
    Fig 1: four subplots against time: CMB temp, fractional inner core radius, magnetic Reynolds number, rate of inner core growth
    Fig 2: absolute value of d^2r/dt^2
    Fig 3: Cr - dr/dt/dTc/dt ith time
            
"""
#import modules
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt 

#Constants
G=6.67e-11 # gravitational constant [N kg^-2 m^2]
year=60*60*24*365 #number of seconds in a year
Myr=1e6*year

#Parameters (taken from Nimmo 2009 - N09)
rho=7019 #[kg m^-3]
cp=835 #[J kg^-1 K^-1]
k=30 #[Wm^-1 K^-1]
alpha=9.2e-5 #[K^-1]
drho=0.05 # \delta rho/rho 
rc=1e5# core radius [m] 
Delta=1.2 #dTm/dP/dT/dP 

#Parameters (my choices)
#change to arrays later
f0=0.01 #initial fractional inner core radius [N09 use fixed value of 0.5]
phiv=1e-12# volumetric ohmic dissipation [WK^-1m^-3] [0.2-2 x1e-12 from N09]
Tc0=np.array([2000]) # intial CMB temp [K]

#Derived constants
D= np.sqrt(3*cp/(2*np.pi*alpha*rho*G)) # scale height [m]
B=D**2/(2*rc**2*(Delta-1))

#import functions
from F_def import F_calc
from f_def import f_calc

# create dTc/dt=g(Tc,f)
def dTcdt(t,Tc):
    """
    

    Parameters
    ----------
    Tc : float
        CMB temp
    Tc0 : float# run the RK method
#teval=np.linspace(0,1,100)
        initial CMB temp
    f0 : float
        intial inner core fraction.
    B : float
        simplified constant D^2/(2rc^2(Delta-1).

    Returns
    -------
    float
       dTc/dt

    """
  
    f=f_calc(Tc,Tc0,f0,B)
    F=F_calc(f)
    return -Tc**2/F*(Delta-1)*phiv/(3*np.pi*rho**2*drho*G*D**2)

# create dridt
def drdt_calc(Tc,dTcdT,f):
    """
    

    Parameters
    ----------
    Tc : float
        core temperature.
    dcdT : float
        core cooling rate (-ve indicates cooling)
    f : float
        inner core size

    Returns
    -------
    float
        rate of inner core growth, dri/dt

    """
    return -D**2*dTcdT/(2*Tc*f*rc*(Delta-1))

#create zero event to stop calc when f=1
def event(t,Tc):
    f=f_calc(Tc,Tc0,f0,B)

    return f-1

event.terminal=True

# run the RK method
#teval=np.linspace(0,1,100)
t_end=20*Myr
step_m=0.001*t_end
yout=scint.solve_ivp(dTcdt,[0,t_end],Tc0,max_step=step_m,events=event)
Tc_solve=yout['y'][0,:]
t_solve=yout['t']
f_out=f_calc(Tc_solve,Tc0,f0,B)

#Calculate the magnetic Reynolds number at each timestep
from Rem_calc import Rem
dTcdt_plot=dTcdt(t_solve,Tc_solve)
drdt=drdt_calc(Tc_solve,dTcdt_plot,f_out)
Re=Rem(rc,drdt)
drdt_plot=drdt*Myr/1e3 #convert to km/Myr from m/s

plt.figure(tight_layout=True)
plt.subplot(2,2,1)
plt.plot(yout['t']/Myr,yout['y'][0,:])
plt.ylabel('CMB T /K')
plt.xlabel('t/Myr')
plt.subplot(2,2,2)
plt.plot(yout['t']/Myr,f_out)
plt.xlabel('t/Myr')
plt.ylabel('f')
plt.subplot(2,2,3)
plt.plot(yout['t']/Myr,Re)
plt.xlabel('t/Myr')
plt.ylabel('Rem')
plt.hlines(y=40,xmin=0,xmax=max(t_solve)/Myr,color='k',linestyle='--')
plt.subplot(2,2,4)
plt.plot(yout['t']/Myr,drdt_plot)
plt.xlabel('t/Myr')
plt.ylabel('dri/dt /km/Myr')
plt.yscale('log')
#plt.savefig('Plots/subplots.png')

#Analytic solution for drdt
def dfdt(f):
    A=phiv/(6*np.pi*rho**2*drho*G*rc**2)
    return A*Tc0/(f*F_calc(f))*np.exp((f0**2-f**2)/(2*B))

#Second derivative -done numerically as kept mucking up analytic differentiation
delta_t=np.ediff1d(t_solve)
def d2fdt2(f):
    return (dfdt(f[1:])-dfdt(f[:-1]))/delta_t
    
drdt_an=dfdt(f_out)*rc*Myr/1e3
sdiv=d2fdt2(f_out)

plt.figure()
plt.plot(t_solve[:373]/Myr,-sdiv[:373],label='-ve d2fdt')
plt.plot(t_solve[374:]/Myr,sdiv[373:],label='+ve d2fdt')
plt.xlabel('t/Myr')
plt.ylabel('|$d^2f/dt^2$| /$Myr^{-2}$')
plt.yscale('log')
plt.savefig('Plots/second_derivative.png')
#plt.legend()

plt.figure()
plt.plot(t_solve[:-1]/Myr,-drdt_an[:-1]/dTcdt_plot[:-1])
plt.xlabel('t/Myr')
plt.ylabel('$C_r$')
plt.yscale('log')
plt.text(8,1e17,'$C_r$=$\\frac{dr}{dt}$/$\\frac{dTc}{dt}$')
plt.savefig('Plots/cr.png')