#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField_f, double *collideField_g, unsigned int* flagField, const double * const wallVelocity,dimensions dim);

void get_sum_of_weights(int x, int y, int z, dimensions d, unsigned int* flagField, double* weight_sum);

void get_Boundary_Temperature(int x, int y, int z, double* T_d);

#endif