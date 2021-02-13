#!/usr/bin/env python3

import numpy as np
import parmed as pmd 
def ambertop_to_grotop(mytop,myincprd,newfilename,mygro):
    parm = pmd.load_file(mytop,myincprd)
    parm.save(newfilename, format='gromacs')
    parm.save(mygro)

def Temperature_calc_exp(T0, k, rep):
    '''
    Calculator for list of temperatures distributed exponentially for use in REMD simulations.
    T0 = Temperature of interest
    k = constant for determining desired acceptance rate* 
    rep = number of replicas
    
    *See Alexandra Patriksson and David van der Spoel, 
    A temperature predictor for parallel tempering simulations, 
    Phys. Chem. Chem. Phys., 10 pp. 2073-2077 (2008) http://dx.doi.org/10.1039/b716554d.
    
    Example
    -----------
    from MDConsole import sysprep as sp
    T0 = 300 # our temperature of interest
    N_rep = 16 # Number of Replicas
    k = 0.005 # constant for determining acceptance rate based on DOF, 
    # value k = 0.005 gives acceptance rate of ~ 20-28% 
    # for a system with 50,000-100,000 atoms simulated without restraints.
    
    >>> myT = anl.Temperature_calc_exp(T0,k,N_rep)
    [300,
     301.5,
     303.02,
     304.53,
     306.06,
     307.59,
     309.14,
     310.69,
     312.24,
     313.81,
     315.38,
     316.96,
     318.55,
     320.15,
     321.75,
     323.37]
    -----------
    '''
    temp_dist = []
    temp_dist.append(T0)
    for i in range(1,rep,1):
        a = T0 * np.exp(k*i)
        b = np.round(a,2) # up to two decimal points
        temp_dist.append(b)
    return temp_dist
