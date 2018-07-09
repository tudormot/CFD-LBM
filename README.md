========================================================================================================

CFD-LAB
SUMMER TERM 2018
GROUP H

========================================================================================================

**LATTICE BOLTZMANN SIMULATION WITH HEAT TRANSFER**

The Heat Transfer implementation followed the procedure proposed by Chih-Hao et al. in "Thermal boun-
dary conditions for thermal lattice Boltzmann simulations"

========================================================================================================

CODE USAGE

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
        - One can modify its parameters in convection.dat.
        - To run:   ./lbsim convection.dat

OUTPUT FOLDERS

- Each simulation is creating an output folder for .vtk files.
- The name of the folder is the name of the scenario.
- The output folder do NOT have to be cleaned, removed or manually created before any simulation.

========================================================================================================
