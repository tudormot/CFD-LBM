#include "collision.h"
#include "LBDefinitions.h"
#include <stdio.h>

void computePostCollisionDistribution(double *currentCell, const double * const tau, const double *const feq){
  
	int Q = 19; // TODO: again...

	for(int i = 0; i < Q; ++i){
		currentCell[i] = currentCell[i] - (1.0/(*tau))*(currentCell[i] - feq[i]);
	}
}

void doCollision(double *collideField, unsigned int *flagField,const double * const tau,dimensions dim, double* vel){
  
	int Q = NO_OF_LATTICE_DIMENSIONS;
	double density;
	double velocity[3];
	double feq[Q];
	double* currentCell;
	int i = 0; // Data structure. Should be pointing to first distribution function
	int xl = dim.xlen +2;
	int xlyl = xl*(dim.ylen+2);


	for(int z=1; z <= dim.zlen; ++z)
		for(int y=1; y <= dim.ylen; ++y)
			for(int x=1; x <= dim.zlen; ++x){ // x is last for efficiency (cache)


				if( flagField[z*xlyl + y*xl + x] == 0){ // index corrected
					currentCell = &collideField[Q*(z*xlyl + y*xl + x) + i];
					computeDensity(currentCell, &density);
					computeVelocity(currentCell, &density, velocity);

					vel[3*(z*xlyl + y*xl + x) + 0] = velocity[0];
					vel[3*(z*xlyl + y*xl + x) + 1] = velocity[1];
					vel[3*(z*xlyl + y*xl + x) + 2] = velocity[2];


					computeFeq(&density, velocity, feq);
					computePostCollisionDistribution(currentCell, tau, feq);
				}
		}

}

