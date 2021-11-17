#!/bin/sh
#SBATCH --job-name=gmx
#SBATCH --output=gpu%j.o
#SBATCH --error=gpu%j.e
#SBATCH --mail-type=END
#SBATCH --mail-user=db4271@nyu.edu
#SBATCH --time=10-00:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=50GB
export MODULEPATH=$MODULEPATH:/xspace/sungroup/modules/
module load gromacs/2020.4_avx2

# Set some environment variables 
SOLVENT=methanol
DIR=/xspace/db4271/dipole_flip/SolventBox
MAINDIR=$DIR/$SOLVENT
echo "Main directory set to $MAINDIR"
MDP=$DIR/MDP
echo ".mdp files are stored in $MDP"

# Change to the location of your GROMACS-2018 installation
GMX=/xspace/sungroup/software/gromacs2020.4_avx2/bin
cd $MAINDIR


    #######################
    # Energy Minimization #
    #######################
echo "Starting energy minimization for solvent = $SOLVENT ..."


$GMX/gmx grompp -f em.mdp -c methanol.gro -p methanol.top -o em_methanol.tpr

$GMX/gmx mdrun -deffnm em_$SOLVENT


    sleep 10
   #######################
    # NVT Equilibration   #
    #######################
echo "Starting constant volume equilibration for solvent = $SOLVENT ..."


$GMX/gmx grompp -f nvt.mdp -c em_methanol.gro -p methanol.top -o nvt_methanol.tpr

$GMX/gmx mdrun -deffnm nvt_$SOLVENT

    #####################
    # NPT EQUILIBRATION #
    #####################

    #################
    # PRODUCTION MD #
    #################

    # End
    echo "Ending. Job completed for solvent = $SOLVENT"

exit;
