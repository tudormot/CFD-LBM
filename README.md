CFD-LAB
SUMMER TERM 2018
GROUP H

**LATTICE BOLTZMANN SIMULATION WITH HEAT TRANSFER**

Simulation of two scenarios following the Lattice Boltzmann Method, taking the BGK approximation for the collision operator.

- Building the code with the Makefile creates the lbsim executable.
- There are two possible scenarios:
    - Driven Cavity: 
        - Classic cavity scenario in 3D, with one moving wall in the upper z-face, 
          free-slip conditions in the y-faces and no-slip conditions on the x-faces 
          and the bottom z-face. 
        - One can modify its parameters in cavity.dat.
        - To run:   ./lbsim cavity.dat
        
    - Natural Convection: 
        - Also a 3D scenario. The velocity boundary conditions are free-slip in the y-faces 
          and no-slip on the rest. The temperature boundary conditions are Dirchlet to the 
          x-faces and adiabatic for the rest.
        - To run:   ./lbsim convection.dat


2.2. Output folders
- Each simulation is creating and output folder for .vtk files.
- Output folder's name consists of chosen scenario and content of "problem" variable.
- Output folders do NOT have to be cleaned, removed or manually created before any simulation.