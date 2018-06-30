#include "collision.h"
#include <stdio.h>

void computePostCollisionDistribution(double *currentCell, const double * const tau, const double *const feq){
  
	int Q = 19; // TODO: again...

	for(int i = 0; i < Q; ++i){
		currentCell[i] = currentCell[i] - (1.0/(*tau))*(currentCell[i] - feq[i]);
	}
}

void doCollision(double *collideField, unsigned int *flagField,const double * const tau,int xlength, double* vel){
  
	int Q = 19;	 // TODO: Q is hardcoded;
	double density;
	double velocity[3];
	double feq[Q];
	double* currentCell;
	int i = 0; // TODO: is this 0 or 1? 5. Data structure. Should be pointing to first distribution function
	int xl2 = (xlength+2)*(xlength+2); // index corrected (worksheet was wrong (?))


	for(int z=1; z <= xlength; ++z)
		for(int y=1; y <= xlength; ++y)
			for(int x=1; x <= xlength; ++x){ // x is last for efficiency (cache)


				if( flagField[z*xl2 + y*(xlength+2) + x] == 0){ // index corrected
					currentCell = &collideField[Q*(z*xl2 + y*(xlength+2) + x) + i];
					computeDensity(currentCell, &density);
					computeVelocity(currentCell, &density, velocity);

					vel[3*(z*xl2 + y*(xlength+2) + x) + 0] = velocity[0];
					vel[3*(z*xl2 + y*(xlength+2) + x) + 1] = velocity[1];
					vel[3*(z*xl2 + y*(xlength+2) + x) + 2] = velocity[2];


					computeFeq(&density, velocity, feq);
					computePostCollisionDistribution(currentCell, tau, feq);
				}
		}

}

