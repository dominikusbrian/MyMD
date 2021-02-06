#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --time=1-12:00:00
#SBATCH --mem=5gb
#SBATCH --job-name="Get-DATA"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=parallel
##SBATCH --qos=argon
##SBATCH --constraint=9242

CASE=1100
min=$((CASE+1))
max=$((CASE+2))

SCRIPTDIR=/xspace2/db4271/Triad/MD_TRAJ/src
TRAJ_STATE="GR"
ENSEMBLE="NVE"
for (( i = $min ; i <= $max ; i++ ))
do
	for REC_STATE in CT1 CT2 GR PI
	do
    REC_FILE="triad_thf_REC_${REC_STATE}_${TRAJ_STATE}_${ENSEMBLE}_$i.dat"
    E_FILE="triad_thf_${REC_STATE}_TRAJ_${TRAJ_STATE}_${ENSEMBLE}_$i.dat"
    awk '{print $2}' $REC_FILE | sed 1d  > $E_FILE
    OUT_FILE="triad_thf_${TRAJ_STATE}_TRAJ_${ENSEMBLE}_$i.out"
    E_FILE="ALL_triad_thf_TRAJ_${TRAJ_STATE}_${ENSEMBLE}_$i.dat"
    awk -f ${SCRIPTDIR}/get-everything.awk $OUT_FILE > $E_FILE
    OUT_FILE="triad_thf_${TRAJ_STATE}_TRAJ_${ENSEMBLE}_$i.out"
    E_FILE="TEMP_triad_thf_TRAJ_${TRAJ_STATE}_${ENSEMBLE}_$i.dat"
    awk -f ${SCRIPTDIR}/get-temperature.awk $OUT_FILE > $E_FILE

	done
done

echo "All $REC_STATE energies on $TRAJ_STATE trajectories ($min ~ $max) were saved"
exit
