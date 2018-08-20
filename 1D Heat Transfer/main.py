# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:49:10 2018

@author: Ahmad

This program numerically solves the 1D Heat Transfer of a metallic rod under
the Dirichlet Boundary Conditions
"""

def Dirichlet():
    # Physical Properties
    kappa = 1 # Thermal conductivity
    rho = 1 # Density of the metal
    Cp = 1 # Heat capacity of the metal
    
    L = 1 # [m]
    
    alpha = kappa/(rho * Cp)
    
    # Boundary and initial conditions
    Ti = 273 # [k]; initial temperature
    T0 = 273 # [k]; boundary temperatures
    
    # Dimensions for discretizations
    Nx = 10