#!/usr/bin/env python3


import numpy as np
import pandas as pd
import ArtisanKit as ak

# For analysis functions

def calc_acceptancerate(filename, whichtemp,print_result):
    '''
    filename = remd log file
    whichtemp = which temperature replica , take scalar 1.0, 2.0 , and so on
    print_result = True / False
    
    Example
    -----------
    >>> path ='/xspace/db4271/MDConsole/MDConsole/REMD'
    >>> f = '%s/rem.log' %(path_rem)
    >>> myaveaccept = calc_acceptancerate(f,1.0,print_result=False)
    >>> calc_acceptancerate(f,1.0,print_result=True)
    
    My acceptance rate =  2.70 %
    
    >>> print(myaveaccept)
    2.7009999999999876
    -----------
    '''
    log= np.loadtxt(filename)
    df = pd.DataFrame(log)
    temp= df.loc[df[0] == whichtemp]
    #rem.log column format:
    #Column 2 = current temperature
    #Column 3 = potential energy
    #Column 6 = success rate

    ave_acceptance = np.mean(temp[6]) * 100
    temp_traj = temp[2]
    Epot = temp[3]
    
    # If zero is ignored
    
    #temp[temp == 0] = np.nan
    #print(temp[6])
    
    #ave_acceptance = np.nanmean(temp.loc[6]) * 100
    if(print_result == True):
        ak.iostream.print_dec("My acceptance rate = ",ave_acceptance, "%", 2)
    else:
        return ave_acceptance , temp_traj, Epot
