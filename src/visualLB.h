#ifndef _VISUALLB_H_
#define _VISUALLB_H_
#include <stdio.h>
#include "LBDefinitions.h"

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(const double* const Vels, const double* const Temps, const unsigned int * const flagField, const char* filename, int t, dimensions dim);
void write_vtkHeader( FILE *fp, dimensions dim);
void write_vtkPointCoordinates( FILE *fp, dimensions dim);

#endif
