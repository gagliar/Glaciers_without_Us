Run a transient simulation with a constant climate (look at a steady solution)

As this simulation restart from 1\_IMPORT\_DEM, the mesh should be the same and the importation of the mesh already executed 

Compile the Surface Mass Balance solver:<br>
`elmerf90 ../SRC/TransientMassBalance.f90 -o bin/MassBalance`<br>

Execute the simulation:<br>
`oarsub -S ./MyRun.sh`

Take care that the numper of nodes required is compatible with the number of partitions of the mesh 
