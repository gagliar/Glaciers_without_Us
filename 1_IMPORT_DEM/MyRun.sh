#!/bin/bash
#################################
#OAR --name 1_INIT_DEM
#OAR --project elmerice
### request ressources: adapt to your needs
#OAR -l /nodes=1,walltime=00:10:00
## request only nodes with 32 cores
#OAR -p n_cores=4 
#################################
#
# Exit on error
set -e

#################################
# Calcul du nombre de noeuds, et nombre de cores
nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`

echo "############"
echo "Nbre de noeuds: " $nbnodes
echo "Nbre de coeurs: " $nbcores
echo "############"
echo 
#################################


### let's go
ulimit -s unlimited
export OMP_NUM_THREADS=1

mpirun -np $nbcores ElmerSolver_mpi

