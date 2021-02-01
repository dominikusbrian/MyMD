#!/usr/bin/env python3

def ambertop_to_grotop(mytop,myincprd,newfilename,mygro):
    parm = pmd.load_file(mytop,myincprd)
    parm.save(newfilename, format='gromacs')
    parm.save(mygro)
