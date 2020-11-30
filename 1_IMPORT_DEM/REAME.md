In this first step, only read the different DEMs needed for the simulation

You can change the resolution of the mesh directly by modifying the **lc** variable in Mesh2d.geo

Then make the Elmer mesh:<br>
`gmsh -1 -2 Mesh2d.geo`<br>
`ElmerGrid 14 2 Mesh2d.msh -autoclean -metis 4 4`<br> 

Execute the simulation:<br>
`./MyRun.sh`

Take care that the numper of nodes required is compatible with the number of partitions of the mesh 
