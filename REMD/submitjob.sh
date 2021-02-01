#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 28
#SBATCH --time=2-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name="test"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=parallel
##SBATCH --qos=argon

source /etc/profile.d/modules.sh
module purge
module load amber
export AMBERHOME=/gpfsnyu/packages/amber/20icc

mpirun -np 28 sander.MPI -O -i 01_min.in -o min.out -p triad_thf_bent_GR.prmtop -c triad_thf_bent.inpcrd -r min.rst
