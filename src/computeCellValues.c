#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
  
	int Q = 19; // TODO: again hardcoded
	(*density) = 0; // reset density, before summing for each cell

	for(int i = 0; i < Q; ++i){
		(*density) += currentCell[i];
	}
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
  
  	double product[3] = {0,0,0}; // reset f*c_i product before new cell calculation
	int Q = 19;
	for(int i = 0; i < Q; ++i){
		product[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
		product[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
		product[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
	}

	velocity[0] = product[0]/(*density); // density still holds sum from previous function
	velocity[1] = product[1]/(*density);
	velocity[2] = product[2]/(*density);


}

void computeFeq(const double * const density, const double * const velocity, double *feq){
 
	int Q = 19; //TODO: hardcoded again -> move it to LBDefinitions(?)
	double cu; // dot product of lattice velocity and velocity
	double cu2; // cu squared;
	double uu; // dot product of velocity with itself
	double cs2 = C_S * C_S; // to avoid squaring in for loop
	double cs4 = cs2 * cs2; // C_S to the power of 4


	for(int i = 0; i < Q; ++i){
		cu = LATTICEVELOCITIES[i][0]*velocity[0] + LATTICEVELOCITIES[i][1]*velocity[1] + LATTICEVELOCITIES[i][2]*velocity[2];
		cu2 = cu*cu;
		uu = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];

		feq[i] = LATTICEWEIGHTS[i]*(*density)*(1 + (cu/cs2) + (cu2/(2*cs4)) - (uu/(2*cs2)));
	}
}

