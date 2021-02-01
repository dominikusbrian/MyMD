#!/bin/bash
#SBATCH -N 4
#SBATCH -n 96
##SBATCH --time=2-00:00:00
#SBATCH --ntasks-per-node=24
#SBATCH --mem=100gb
#SBATCH --job-name="eq-REMD"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=argon
#SBATCH --qos=argon
#SBATCH --constraint=g6146

source /etc/profile.d/modules.sh
module purge
module load amber/20z

mpirun -np 96 pmemd.MPI -ng 8 -groupfile equilibrate.groupfile 
