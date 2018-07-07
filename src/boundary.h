#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField_f, double *collideField_g, unsigned int* flagField, const double * const wallVelocity,dimensions dim, double T_cold, double T_warm);

void get_sum_of_weights(int x, int y, int z, dimensions d, unsigned int* flagField, double* weight_sum);

void get_Boundary_Temperature(int flag, double* T_d, double T_cold, double T_warm);

int look_up_vel(int dx, int dy, int dz);

#endif
