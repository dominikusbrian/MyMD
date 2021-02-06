
#Global Parameters
TRAJ_STATE="NVT"
Donor="CT1"

for CASE in {1100..22000..100}
do
for CONFIG in FB
do

SRCDIR=/xspace2/db4271/Triad/MD_TRAJ/src
JOBDIR=/xspace/db4271/Triad/MD_TRAJ/${CONFIG}/${Donor}-${TRAJ_STATE}-5fs/${Donor}-${TRAJ_STATE}-${CASE}
RSTDIR=/xspace2/db4271/Triad/MD_TRAJ/${CONFIG}/HOMEBASE/SAMPLING_${Donor}_NVT
mkdir -p ${JOBDIR}

#Take RST for the production
cd ${RSTDIR}
cp *_${CASE}.rst ${JOBDIR}

#Setup slurm script for the run
cd ${SRCDIR}
cp amber20z_gpu-runlong.slurm  ${JOBDIR} 
cd ${JOBDIR}
sed_param=s/CASE=.*/CASE="${CASE}"/
sed -i "$sed_param" amber20z_gpu-runlong.slurm
sed_param=s/CONFIG=.*/CONFIG="${CONFIG}"/
sed -i "$sed_param" amber20z_gpu-runlong.slurm

echo "slurm script successfully edited, submitting now"

sbatch amber20z_gpu-runlong.slurm

done
done

