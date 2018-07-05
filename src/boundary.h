#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField_f,double *collideField_g, unsigned int* flagField, const double * const wallVelocity,dimensions dim);

#endif

