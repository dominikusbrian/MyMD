#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 2
#SBATCH --time=1-12:00:00
#SBATCH --mem=50gb
#SBATCH --job-name="Get-DATA"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=parallel
## TOTAL = 200 Trajectories ##
TRAJ_STATE="NVT"
DONOR="GR" 

for i in {1101..21001..100} #200 TRAJ
do 
	#cat TEMP_triad_thf_TRAJ_${DONOR}_${TRAJ_STATE}_${i}.dat >> TEMP.dat
	for REC_STATE in GR PI CT1 CT2
	do
        cat triad_thf_${REC_STATE}_TRAJ_${DONOR}_${TRAJ_STATE}_${i}.dat >> TRIAD_E_${REC_STATE}_TRAJ_${DONOR}_${TRAJ_STATE}.dat
        done
	
done
