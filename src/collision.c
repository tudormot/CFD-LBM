#include "collision.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include "helper.h"
#include <math.h>

void computePostCollisionDistribution(double *currentCell_f, double *currentCell_g, const double * const tau_f, const double * const tau_g,
		const double *const feq, const double * const geq, double* F_b){
  
	int Q = NO_OF_LATTICE_DIMENSIONS;

	for(int i = 0; i < Q; ++i){
		currentCell_f[i] = currentCell_f[i] - (1.0/(*tau_f))*(currentCell_f[i] - feq[i]);// + F_b[i]; // Gravity is probably missing here
		currentCell_g[i] = currentCell_g[i] - (1.0/(*tau_g))*(currentCell_g[i] - geq[i]); //
	}
}

void doCollision(double *collideField_f, double *collideField_g, unsigned int *flagField,const double * const tau_f, const double * const tau_g ,
		dimensions dim, double* vel, double* Temps){
  
	int Q = NO_OF_LATTICE_DIMENSIONS;
	double density;
	double velocity[3];
	double feq[Q];
	double geq[Q];
	double* currentCell_f;
	double* currentCell_g;
    double Temp;
    double F_b[Q];
	int xl = dim.xlen +2;
	int xlyl = xl*(dim.ylen+2);


	for(int z=1; z <= dim.zlen; ++z)
		for(int y=1; y <= dim.ylen; ++y)
			for(int x=1; x <= dim.zlen; ++x){ // x is last for efficiency (cache)


				if( is_fluid(flagField[z*xlyl + y*xl + x])){
					currentCell_f = &collideField_f[Q*(z*xlyl + y*xl + x)];
					computeDensity(currentCell_f, &density);
					computeVelocity(currentCell_f, &density, velocity);

					vel[3*(z*xlyl + y*xl + x) + 0] = velocity[0];
					vel[3*(z*xlyl + y*xl + x) + 1] = velocity[1];
					vel[3*(z*xlyl + y*xl + x) + 2] = velocity[2];
					currentCell_g = &collideField_g[Q*(z*xlyl + y*xl + x)];
					computeTemperature(currentCell_g, &density, &Temp);
					Temps[z*xlyl + y*xl + x] = Temp; //check index

					computeFeq(&density, velocity, feq);
					computeGeq(feq, &Temp, geq);
					computeBuoyancy(&Temp, F_b);
					computePostCollisionDistribution(currentCell_f, currentCell_g, tau_f, tau_g, feq, geq, F_b);
				}
		}

}

void computeBuoyancy(double* Temp, double* F_b){

	double G;
	double beta = 0.00021;
	double T_c = 0;
	double T_w = 10;
	double T_m = (T_c + T_w)/2;
	double g = -9.81;
	G = beta*g*(*Temp - T_m);

	for(int i = 0; i < NO_OF_LATTICE_DIMENSIONS; ++i){
		F_b[i] = 3*LATTICEWEIGHTS[i]*G*LATTICEVELOCITIES[i][2];
	}
}
