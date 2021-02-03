#!/usr/bin/env python3


import numpy as np
import pandas as pd
import ArtisanKit as ak

# For analysis functions

def calc_acceptancerate(filename, whichtemp,print_result):
    '''
	Function to calculate for Acceptance rate from AMBER rem.log file

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
    #Column 6 is the from typical amber format for rem.log file 
    ave_acceptance = np.mean(temp[6]) * 100

    if(print_result == True):
        ak.iostream.print_dec("My acceptance rate = ",ave_acceptance, "%", 2)
    else:
        return ave_acceptance
