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
SOLVENT=THF
DIR=/xspace/db4271/dipole_flip/SolventBox
MAINDIR=$DIR/$SOLVENT/result
echo "Main directory set to $MAINDIR"
MDP=$DIR/MDP
echo ".mdp files are stored in $MDP"

# Change to the location of your GROMACS-2018 installation
GMX=/xspace/sungroup/software/gromacs2020.4_avx2/bin
cd $MAINDIR

    ##############################
    # ENERGY MINIMIZATION STEEP  #
    ##############################
    echo "Starting minimization for solvent = $SOLVENT..." 

    mkdir EM
    cd EM

    # Iterative calls to grompp and mdrun to run the simulations

    $GMX/gmx grompp -f $MDP/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/methane_water.gro -p $FREE_ENERGY/topol.top -o min$LAMBDA.tpr

    $GMX/gmx mdrun -deffnm min$LAMBDA

    sleep 10

    #####################
    # NVT EQUILIBRATION #
    #####################



    #####################
    # NPT EQUILIBRATION #
    #####################




    #################
    # PRODUCTION MD #
    #################



    # End
    echo "Ending. Job completed for solvent = $SOLVENT"
exit;
