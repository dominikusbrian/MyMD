#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 24
#SBATCH --time=2-00:00:00
#SBATCH --mem=5gb
#SBATCH --job-name="file cut "
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=argon
#SBATCH --qos=argon


CASE=2
for (( j = 1; j <= 5; j++))
do
blocksize=8000000
tblocksize=80
beforej=$((j-1))
BASIS=$((beforej*blocksize))
START=$((1+BASIS ))
STOP=$((j*blocksize))
TBASIS=$((beforej*tblocksize))
TSTART=$((1+TBASIS ))
TSTOP=$((j*tblocksize))

if [ $CASE  == ${j} ]
	then
	sed -i ${TSTART},${TSTOP}!d TEMP.dat
	for STATE in GR PI CT1 CT2
	do   		
	sed -i ${START},${STOP}!d TRIAD_E_${STATE}_TRAJ_PI_NVE.dat
	done
fi	
done
