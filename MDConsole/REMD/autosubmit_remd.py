#!/usr/bin/env python3

import os 
os.system("echo Hello from the other side!")

#setup_equil()


def slurm_header(partition, N_nodes, N_cores, restriction1, restriction2, restriction3):
    # "#!/bin/bash
    # "#SBATCH -N N_nodes
    # "#SBATCH -n 96
    # "##SBATCH --time=2-00:00:00
    # "#SBATCH --ntasks-per-node=24
    # "#SBATCH --mem=100gb
    # "#SBATCH --job-name="REMD"
    # "#SBATCH --mail-type=END
    # "#SBATCH --mail-user=db4271@nyu.edu
    # "#SBATCH --output=%job.o
    # "#SBATCH --error=%job.e
    # "#SBATCH --partition=argon
    # "#SBATCH --qos=argon
    # "#SBATCH --constraint=g6146
    pass


def run_energy_minimization(self):
    """
    Running energy minimization.
    Take min.in as input file and produce a set of outputs:
    min.out
    min.rst 
    """
    print(self)
    os.system("echo submitting energy minimization job ")
    os.system("sbatch submitjob.sh")

def setup_equil():
    os.system("echo executing setup_equilibrate_input.x")
    os.system("./setup_equilibrate_input.x")
    return 

def run_heat_and_equil():
    os.system("echo submitting heating and equilibrium job ")
    os.system("sbatch submitjob_mpi.sh")

def setup_remd():
    os.system("echo executing setup_remd_input.x")
    os.system("./setup_remd_input.x")
    return 

def run_remd():
    os.system("echo submitting remd simulation job ")
    os.system("sbatch submitjob_remd.sh")

