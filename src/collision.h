#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"
#include "LBDefinitions.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistribution(double *currentCell_f, double *currentCell_g, const double * const tau_f, const double * const tau_g,
		const double *const feq, const double * const geq, double* F_b);


/** carries out the whole local collision process. Computes density, velocity, temperature and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(double *collideField_f, double *collideField_g, unsigned int *flagField,const double * const tau_f, const double * const tau_g, dimensions dim, double* vel, double* Temps,
		double beta, double  T_cold, double T_warm, double gravity);

/*computes F_b, term which models contribution of temperature change + gravity to the density field*/
void computeBuoyancy(double* Temp, double* F_b, double beta, double T_cold, double T_warm, double gravity);
#endif

