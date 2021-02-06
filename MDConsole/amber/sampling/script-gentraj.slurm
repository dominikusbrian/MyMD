#!/bin/bash
#SBATCH -N 2
#SBATCH -n 21
#SBATCH --time=5-00:00:00
#SBATCH --mem=30gb
#SBATCH --job-name="CT1-FB-gentraj"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%job.o
#SBATCH --error=%job.e
#SBATCH --partition=parallel
##SBATCH --qos=parallel
##SBATCH --constraint=g6146

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


CONFIG="FB"
TRAJ_STATE="NVT"
Donor="CT2"
Initial_STATE="triad_thf_${Donor}_TRAJ_${TRAJ_STATE}"
WORKDIR=/xspace2/db4271/Triad/${CONFIG}/HOMEBASE/SAMPLING_${Donor}_${TRAJ_STATE}
HOMEDIR=/xspace2/db4271/Triad/${CONFIG}/HOMEBASE
HQDIR=/xspace2/db4271/Triad/src

for j in {21300..22000..100}
do
### recalculate potential energies with input trajectory (netcdf)
# [1] MD simulation
beforej=$((j-100))
mpirun -np 21 pmemd.MPI -O -i ${HQDIR}/NFF-NVT-gentraj-FB.mdin \
        -o ${WORKDIR}/triad_thf_${Donor}_TRAJ_${TRAJ_STATE}_${j}.out \
        -p ${HQDIR}/triad_thf_${CONFIG}_${Donor}.prmtop \
        -c ${WORKDIR}/${Initial_STATE}_${beforej}.rst \
        -ref ${WORKDIR}/${Initial_STATE}_${beforej}.rst \
        -r ${WORKDIR}/triad_thf_${Donor}_TRAJ_${TRAJ_STATE}_${j}.rst \
        -inf ${WORKDIR}/triad_thf_${Donor}_TRAJ_${TRAJ_STATE}_${j}.mdinfo \
        -x ${WORKDIR}/triad_thf_${Donor}_TRAJ_${TRAJ_STATE}_${j}.mdcrd.nc
    echo traj_${j} finished on `date`
done
