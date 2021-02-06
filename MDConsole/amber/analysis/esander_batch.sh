#!/bin/sh
#SBATCH --job-name=GR_REC
##SBATCH --time=1-00:00:00
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH -p argon
#SBATCH --nodes=1
#SBATCH -c 28
#SBATCH --mem=50GB
#SBATCH --qos=argon
#SBATCH --constraint=g6132
##SBATCH --constraint=g6146|9242

# load amber20z module automatically according to cpu model
module purge
export MODULEPATH=$MODULEPATH:/xspace/sungroup/modules
# for e5 series cpu node (compute1-27), load avx2 version
if [ -n ${SLURM_JOB_NODELIST:7} ] && [ ${SLURM_JOB_NODELIST:7} -lt 28 ];then
module load amber/20z_avx2
else
# for gold/platinum seris cpu node, load avx512 version
module load amber/20z
fi


CASE=2600
START=$((CASE+1))
STOP=$((CASE+2))
CONFIG="FB"
TRAJ_STATE="NVT"
DONOR="GR"


JOBDIR=/xspace2/db4271/Triad/${CONFIG}/${DONOR}-${TRAJ_STATE}-5fs/${DONOR}-${TRAJ_STATE}-${CASE}
HQDIR=/xspace2/db4271/Triad/${CONFIG}/scriptHQ


#Energy Calculation with esander
for (( j = ${START}; j <= ${STOP}; j++))
do
for STATE in GR PI CT1 CT2 
do
cat > esander_${STATE}.in <<EOF
parm ${HQDIR}/triad_thf_${CONFIG}_${STATE}.prmtop
trajin triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}_${j}.mdcrd.nc
esander ${DONOR} out triad_thf_REC_${STATE}_${DONOR}_${TRAJ_STATE}_${j}.dat ntb 1 cut 12 ntf 2 ntc 2
EOF
mpirun -np 28 cpptraj.MPI -i esander_${STATE}.in > esander_${STATE}.out
done
done


