# FluidFreeSurface
A Demo for Fluid Simulation with Free Surface

## Features
 * Based on the Marker-And-Cell (MAC) method
 * Used staggered MAC grid
 * Free surface presented
 * Supported both 2D and 3D
 * Used Blobby and Marching Cubes for creating surface mesh
 * Used freeglut (alternative to OpenGL) for real-time rendering
 * Used Raytracing in Maya 2015 for non real-time rendering
 * Programed using C++
 
## Algorithm Framework

For each iteration:

1. Calculate the simulation time step   
2. Update the particles¡¯ position according to the current velocity field  
3. Update each cell type as fluid or air according to whether it contains particles  
4. Set boundary conditions for new fluid cells     
5. Update current velocity field based on the Navier-Strokes Equations
 * Apply advection using "semi- Lagrangian" method
 * Apply external forces
 * Apply viscosity using Gauss-Seidel Relaxtion (diffusion)
 * Apply the press projection by solving a liner system
 * Extrapolate fluid velocities into nearby air cells
 * Set boundary conditions for solid cells and cells on the free surface
  
## References
[Blinn J., "A Generalization of Algebraic Surface Drawing," ACM Trans. Graph., 1982.]
(https://www.microsoft.com/en-us/research/wp-content/uploads/1982/07/p235-blinn.pdf)

[Welch J. E., Harlow F. H., Shannon J.
P., and Daly B. J., "THE MAC METHOD a computing
technique for solving viscous, incompressible, transient
fluid-flow problems involving free surfaces," Report
LA-3425, Los Alamos Scientific Laboratory, 1965.]
(http://www.cs.rpi.edu/~cutler/classes/advancedgraphics/S13/papers/harlow_welch.pdf)

[Mark C., Greg T. and Peter M., "Rigid, Melting, and Flowing Fluid," Doctoral Dissertation, 2004.]
(http://www.eng.utah.edu/~cs6660/carlsonthesis.pdf)

[Bridson R., "Fluid Simulation for Computer Graphics", CRC Press, 2008.]
(https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)

## Links
[Detailed description on my homepage](http://zhanghaotian1994.com/projects/FluidFreeSurface/)
