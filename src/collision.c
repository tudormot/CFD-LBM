#include "collision.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
  /* TODO */
	int Q = 19; // TODO: again...

	for(int i = 0; i < Q; ++i){
		currentCell[i] = currentCell[i] - (1/(*tau))*(currentCell[i] - feq[i])
	}
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
  /* TODO */
	int Q = 19;	 // TODO: Q is hardcoded;
	double density;
	double velocity[3];
	double feq[Q];
	double* currentCell = NULL;
	int i = 0; // TODO: is this 0 or 1? 5. Data structure. Should be pointing to first distribution function
	int xl2 = (xlength+2)*(xlength+2); // index corrected (worksheet was wrong (?))


	for(int z=1; z <= xlength; ++z)
		for(int y=1; y <= xlength; ++y)
			for(int x=1; x <= xlength; ++x){ // x is last for efficiency (cache)


				if( flagField[z*xl2 + y*(xlength+2) + x] == 0){ // index corrected
					currentCell = &collideField[Q*(z*xl2 + y*(xlength+2) + x) + i];
					computeDenisty(currentCell, &density);
					computeVelocity(currentCell, &density, velocity);
					computeFeq(density, velocity, feq);
					computePostCollisionDistribution(currentCell, tau, feq);
				}
		}

}

