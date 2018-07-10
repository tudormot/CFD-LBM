#ifndef _STREAMING_H_
#define _STREAMING_H_
#include "LBDefinitions.h"

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
/* streaming step as detailed in worksheet 2. The function moves cell values from source cell to destination cell,
 * based on the lattice velocities directions*/
void doStreaming(double *collideField_f, double *streamField_f,double * collideField_g, double * streamField_g, unsigned int *flagField,dimensions dim);

#endif

