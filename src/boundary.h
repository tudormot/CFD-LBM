#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"

/** handles the boundaries in our simulation setup */
/* modifies values of boundary cells for both density and energy field.
 * Theory used based on old CFD-worksheet2(for densityfield, both Neumann and Dirichlet) and
 * Chih-Hao et al. in "Thermal boundary conditions for thermal lattice Boltzmann simulations"
 * for energy field(both Neumann and Dirichlet).
 * Function should be called at each timestep*/
void treatBoundary(double *collideField_f, double *collideField_g, unsigned int* flagField, double *Temps, const double * const wallVelocity,dimensions dim, double T_cold, double T_warm);

/*helper function, used in setting boundary conditions for the energy field*/
void get_sum_of_weights(int x, int y, int z, dimensions d, unsigned int* flagField, double* weight_sum);

/*helper function, used in setting boundary conditions for the energy field*/
void get_Boundary_Temperature(int flag, double* T_d, double T_cold, double T_warm);

/*helper function, used in setting boundary conditions for the energy field*/
int look_up_vel(int dx, int dy, int dz);

#endif
