#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
    
    /*declarations*/
    double *collideField=NULL;
    double *streamField=NULL;
    unsigned int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    unsigned int timesteps;
    unsigned int timestepsPerPlotting;
    double* Vels, Temps;
    
    /*read parameters*/
    printf("=================================================================\n");
    printf("\nSIMULATION PARAMETERS:\n\n");
    readParameters(&xlength, &tau,velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);
    
    /*remove old results*/
    system("rm -rf Output");
    system("mkdir Output");
    printf("\n");
    
    /*memory allocation of required arrays*/
    int size = (xlength+2) * (xlength+2) * (xlength+2);
    collideField_f = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double) );    //Collide field for f
    streamField_f = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double) );    //Stream field for f
    collideField_g = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double) );    //Collide field for g
    streamField_g = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double) );    //Stream field for g
    flagField = (unsigned int*) calloc(size, sizeof(int));
    Vels = (double*) calloc(size*3, sizeof(double));
    Temps = (double*) calloc(size, sizeof(double));
    
    /*initialization of fields*/
    initialiseFields(collideField,streamField,flagField,xlength, velocityWall);
    
    /*now start the calculation: */
    printf("=================================================================\n");
    printf("\nCOMPUTING...\n\n");
    for(int t = 1; t <= timesteps; t++){
        
        /*stream and put the results in collideField vector */
        double *swap=NULL;
        doStreaming(collideField,streamField,flagField,xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        
        /*collide*/
        doCollision(collideField,flagField,&tau,xlength,Vels);
        
        /*apply boundary conditions*/
        treatBoundary(collideField,flagField,velocityWall,xlength);
        
        /*write output and print progress to console*/
        if (t%timestepsPerPlotting==0){
            writeVtkOutput(collideField,flagField,argv[0],t,xlength, Vels);
            printf("t = %4d/%d (%.1f%% completed)\n",t,timesteps, 100.0*t/timesteps);
        }
    }
    
    /*free up memory*/
    free(collideField);
    free(streamField);
    free(flagField);
    free(Vels);
    
    printf("\n=================================================================\n");
    printf("SIMULATION ENDED SUCCESFULLY");
	printf("\n=================================================================\n");
    return 0;
}

#endif

