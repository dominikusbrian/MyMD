#!/bin/sh
#SBATCH --job-name=CT1-NVT-5fs
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --time=10-00:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=50GB

# load amber20z module automatically according to gpu model
module purge
export MODULEPATH=$MODULEPATH:/xspace/sungroup/modules
if [ "$SLURM_JOB_NODELIST" == "gpu7" ];then
# for 2080Ti, load avx512_cuda11.02 version
module load amber/20z
elif [ "$SLURM_JOB_NODELIST" == "gpu6" ];then
# for P100, load avx2_cuda11.02 version
module load amber/20z_avx2
else
# for K10, load avx2_cuda10.2 version
module load amber/20z_avx2_cuda10.2
fi

CASE=1100
CONFIG="FB"
TRAJ_STATE="NVT"
DONOR="CT1"

# make a unique diretory in /tmp/ and copy files needed to run md
JOBDIR=/xspace/db4271/Triad/MD_TRAJ/${CONFIG}/${DONOR}-${TRAJ_STATE}-5fs/${DONOR}-${TRAJ_STATE}-${CASE}
TEMPDIR=/tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
HQDIR=/xspace2/db4271/Triad/MD_TRAJ/src
mkdir -p $TEMPDIR
cd ${JOBDIR}
cp *_${CASE}.rst $TEMPDIR
cd $HQDIR
cp *triad_thf_${CONFIG}_${DONOR}.prmtop* ${CONFIG}-NFF-${TRAJ_STATE}-traj.mdin $TEMPDIR
cd $TEMPDIR

Initial_STATE="triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}"
START=$((CASE+1))
STOP=$((CASE+1))

for (( j = ${START}; j <= ${STOP}; j++))
do
# [1] MD simulation sunning in TEMPDIR
beforej=$((j-1))
    pmemd.cuda -O -i ${CONFIG}-NFF-${TRAJ_STATE}-traj.mdin \
        -o triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}_${j}.out \
        -p triad_thf_${CONFIG}_${DONOR}.prmtop \
        -c ${Initial_STATE}_${beforej}.rst \
        -ref ${Initial_STATE}_${beforej}.rst \
        -r triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}_${j}.rst \
        -inf triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}_${j}.mdinfo \
        -x triad_thf_${DONOR}_TRAJ_${TRAJ_STATE}_${j}.mdcrd.nc
    echo traj_${j} finished on `date`
done

# move mdcrd and output files needed to $JOBDIR
mv *.out* *.mdcrd* *.rst $JOBDIR
rm -rf $TEMPDIR
