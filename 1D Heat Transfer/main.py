# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:49:10 2018

@author: Ahmad

This program numerically solves the 1D Heat Transfer of a metallic rod under
the Dirichlet Boundary Conditions
"""
import numpy as np
from scipy.integrate import odeint

def dydt(y,t,T0,alpha,Nx,x):
    dydt = np.zeros(Nx)
    dydt[1:-1] = alpha * (y[2:] - 2*y[1:-1] + y[:-2]) / (np.diff(x)[1:])**2
    
    return dydt

def AnalyticalSol(x,t,kupp,alpha,T0,Ti,L):
    xx,tt = np.meshgrid(x,t)
    res = np.zeros(xx.shape)
    for k in np.arange(0,kupp+1,):
        res += np.sin((2*k + 1)*np.pi*xx/L)*np.exp(-alpha*(2*k + 1)**2*np.pi**2*tt/L**2)/(2*k + 1)
    res = T0 + 4*(Ti - T0)*res/np.pi
    return res

def Dirichlet():
    # Physical Properties (Aluminum)
    kappa = 273 # [W/m-K] Thermal conductivity
    rho = 2702 # [kg/m^3] Density of the metal
    Cp = 903 # [J/kg-K] Heat capacity of the metal
    
    alpha = kappa/(rho * Cp)
    
    # Boundary and initial conditions
    Ti = 273 # [k]; initial temperature
    T0 = 300 # [k]; boundary temperatures
    
    # Dimensions for discretizations
    L = 10 # [m]
    Nx = 30
    x = np.linspace(0,L,Nx)
    
    # initial condition vector
    y0 = np.ones(Nx) * Ti
    y0[0] = T0; y0[-1] = T0
    
    # solution time
    hr = 1 # [h]
    Nt = hr*3600 # view evolution every second
    timeArr = np.linspace(0,hr*3600,Nt) # [s]
    
    sol = odeint(dydt,y0,timeArr,
                 args=(T0,alpha,Nx,x))
    
    ansol = AnalyticalSol(x,timeArr,1000,alpha,T0,Ti,L)
    return sol,ansol
    
sol,ansol = Dirichlet()
print(sol[-1])
print(ansol[-1])