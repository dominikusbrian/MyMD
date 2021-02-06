
#Global Parameters
TRAJ_STATE="NVT"
Donor="GR"

for CASE in {1000..2000..100}
do

for CONFIG in FB
do

SRCDIR=/xspace2/db4271/Triad/MD_TRAJ/src
JOBDIR=/xspace/db4271/Triad/MD_TRAJ/${CONFIG}/${Donor}-${TRAJ_STATE}-LONG/${Donor}-${TRAJ_STATE}-LONG_${CASE}
RSTDIR=/xspace2/db4271/Triad/MD_TRAJ/${CONFIG}/HOMEBASE/SAMPLING_${Donor}_NVT
mkdir -p ${JOBDIR}

#Take RST for the production
cd ${RSTDIR}
cp triad_thf_GR_TRAJ_NVT_3300.rst ${JOBDIR}
MYRST=triad_thf_GR_TRAJ_NVT_3300.rst

#Setup slurm script for the run
cd ${SRCDIR}
cp amber20z_gpu-runsuperlong.slurm  ${JOBDIR} 
cd ${JOBDIR}
cp ${MYRST} triad_thf_${Donor}_TRAJ_${TRAJ_STATE}_${CASE}.rst

sed_param=s/CASE=.*/CASE="${CASE}"/
sed -i "$sed_param" amber20z_gpu-runsuperlong.slurm
sed_param=s/CONFIG=.*/CONFIG="${CONFIG}"/
sed -i "$sed_param" amber20z_gpu-runsuperlong.slurm

echo "slurm script successfully edited, submitting now"

sbatch amber20z_gpu-runsuperlong.slurm

done
done

