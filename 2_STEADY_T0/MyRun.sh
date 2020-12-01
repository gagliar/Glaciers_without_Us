#!/bin/bash
#################################
#OAR --name 1_INIT_DEM
#OAR --project elmerice

#OAR -O INIT_DEM.%jobid%.o
#OAR -E INIT_DEM.%jobid%.e

#OAR -l nodes=1/core=4,walltime=12:00:00
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


