#!/usr/bin/env python3


import sys
sys.path.append(r'/xspace/db4271/MDConsole/')
sys.path.append(r'/xspace/db4271/ArtisanKit/')
import numpy as np
import pandas as pd
import ArtisanKit as ak
import MDAnalysis as mda
import nglview as nv
import pytraj as pt
from scipy.spatial import distance

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


def calculate_structure (top,traj, atomA, atomB, atomC, atomD):
    '''
    top = input topology file
    traj = input MD trajectory
    atomA, atomB, atomC, atomD = 4 atoms selected as physical features
    e2edist = the end-to-end distance is measured from euclidian distance between 
    atomA and atomD.
    angle_1 = angle between A-B-C
    anlge_2 = angle between A-C-D
    
    Usage example
    -----------
    >>> from MDConsole import analysis as anl
    >>> top = '/xspace/db4271/STCT/noTHF/triad_NOthf_bent_GR.prmtop'
    >>> traj = '/xspace/db4271/STCT/noTHF/remd_012_nothf.nc'
    >>> time_GR, rgyr_GR, e2edist_GR, angle_1_GR, angle_2_GR = anl.calculate_structure(top, traj, atomA, atomB, atomC, atomD)

    -----------
    '''

    u = mda.Universe(top,traj, format='NC')
    conf = u.trajectory[:]
    time = []
    rgyr = []
    e2edist = []
    angle_1 = []
    angle_2 = []
    for ts in conf:
        #extract time
        myt = u.trajectory.time
        #calculate gryation radius
        gr = u.atoms.radius_of_gyration()
        # calculate end to end distance
        myL = u.atoms[[atomA, atomD]]
        atom1 = myL.positions[0]
        atom2 = myL.positions[1]
        myL = distance.euclidean(atom1,atom2)
        # extract pseudo angles
        alpha = u.atoms[[atomA, atomB, atomC]]
	beta = u.atoms[[atomA,atomC,atomD]]
        myangle1 = alpha.angle
	myangle2 = beta.angle
        
        #append all the data to the list
        time.append(myt)
        rgyr.append(gr)
        e2edist.append(myL)
        angle_1.append(myangle1.value())
	angle_2.append(myangle2.value())
    return time, rgyr, e2edist, angle_1, angle_2 
