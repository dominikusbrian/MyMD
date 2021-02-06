#!/bin/bash
# Generate Slurm script and submit jobs  for all cases in the Linear Response ( Nonlinear Response ) IMT calculation LRi



for CASE in NLR 
do

for conformation in FB 
do
WORKDIR=/xspace2/db4271/Triad/LR-IMT/
OUTDIR=${WORKDIR}/MD_Data/${conformation}/${CASE}
donor='PI'
donor_id=0 # for donor = PI
ensemble='NVE'
TRAJ='TRAJ_PI-NVE'
#acceptor_id 1 = CT1 , acceptor_id 2 = CT2;
for acceptor_id in 1 2 
do


Code=${WORKDIR}/src/${CASE}_IMT_${conformation}
energy_path=${WORKDIR}/MD_Data/${conformation}/${TRAJ}

mkdir -p ${OUTDIR}/${donor}CT${acceptor_id}
cd ${energy_path} 
cp *TRIAD* TEMP.dat ${OUTDIR}/${donor}CT${acceptor_id}

cd ${WORKDIR}/src
cp ${CASE}_IMT_${conformation} ${OUTDIR}/${donor}CT${acceptor_id}

cd ${OUTDIR}/${donor}CT${acceptor_id}

cat > submit-${CASE}_IMT_${conformation}_${donor}CT${acceptor_id}.slurm <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --time=1-12:00:00
#SBATCH --mem=5gb
#SBATCH --job-name="IMT"
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --partition=parallel

${Code} ${donor_id} ${acceptor_id}
EOF
sbatch submit-${CASE}_IMT_${conformation}_${donor}CT${acceptor_id}.slurm
done
done
done
