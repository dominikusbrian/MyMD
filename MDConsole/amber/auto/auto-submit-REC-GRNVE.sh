
#Global Parameters
TRAJ_STATE="NVE"
Donor="GR"

for CASE in {20000..20000..100}
do

for CONFIG in FL
do

SRCDIR=/xspace2/db4271/Triad/MD_TRAJ/src
JOBDIR=/xspace2/db4271/Triad/MD_TRAJ/${CONFIG}/${Donor}-${TRAJ_STATE}-5fs/${Donor}-${TRAJ_STATE}-${CASE}
mkdir -p ${JOBDIR}


#Setup slurm script for the run
cd ${SRCDIR}
cp esander_batch.sh  ${JOBDIR} 
cd ${JOBDIR}
sed_param=s/CASE=.*/CASE="${CASE}"/
sed -i "$sed_param" esander_batch.sh
sed_param=s/CONFIG=.*/CONFIG="${CONFIG}"/
sed -i "$sed_param" esander_batch.sh

echo "slurm script successfully edited, submitting now"

sbatch esander_batch.sh

done
done

