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
    echo "Starting constant volume equilibration..."

    cd ../
    mkdir NVT
    cd NVT

    $GMX/gmx grompp -f $MDP/nvt_$LAMBDA.mdp -c ../EM/min$LAMBDA.gro -p $FREE_ENERGY/topol.top -o nvt$LAMBDA.tpr

    $GMX/gmx mdrun -deffnm nvt$LAMBDA

    echo "Constant volume equilibration complete."

    sleep 10

    #####################
    # NPT EQUILIBRATION #
    #####################
    echo "Starting constant pressure equilibration..."

    cd ../
    mkdir NPT
    cd NPT

    $GMX/gmx grompp -f $MDP/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/topol.top -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr

    $GMX/gmx mdrun -deffnm npt$LAMBDA

    echo "Constant pressure equilibration complete."

    sleep 10

    #################
    # PRODUCTION MD #
    #################
    echo "Starting production MD simulation..."

    cd ../
    mkdir Production_MD
    cd Production_MD

    $GMX/gmx grompp -f $MDP/md_$LAMBDA.mdp -c ../NPT/npt$LAMBDA.gro -p $FREE_ENERGY/topol.top -t ../NPT/npt$LAMBDA.cpt -o md$LAMBDA.tpr

    $GMX/gmx mdrun -deffnm md$LAMBDA

    echo "Production MD complete."

    # End
    echo "Ending. Job completed for lambda = $LAMBDA"

    cd $FREE_ENERGY
exit;
