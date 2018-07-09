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
    double *collideField_f=NULL;
    double *streamField_f=NULL;
    double *collideField_g=NULL;
    double *streamField_g=NULL;
    unsigned int *flagField=NULL;
    dimensions dim;    //struct that contains the domain dimensions
    double tau_f, tau_g, T_cold, T_warm, beta, gravity;
    double velocityWall[3];
    unsigned int timesteps;
    unsigned int timestepsPerPlotting;
    double* Vels;
    double* Temps;
    char * filename= NULL; //initialised in readParameters, will point to filename
    
    /*read parameters*/
    printf("=================================================================\n");
    printf("\nSIMULATION PARAMETERS:\n\n");

    readParameters(&dim, &tau_f, &tau_g, &T_cold, &T_warm, &beta, &gravity, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv, &filename);
    printf("\n --------------------->beta: %f,  gravity: %f \n", beta, gravity);
    
    /*remove old results*/
    system("rm -rf Output");
    system("mkdir Output");
    printf("\n");
    
    /*memory allocation of required arrays*/
    int size = (dim.xlen+2) * (dim.ylen+2) * (dim.zlen+2);
    collideField_f = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double));
    streamField_f = (double*) calloc( size * NO_OF_LATTICE_DIMENSIONS, sizeof(double));
    collideField_g = (double*) calloc(size * NO_OF_LATTICE_DIMENSIONS, sizeof(double));
    streamField_g = (double*) calloc( size * NO_OF_LATTICE_DIMENSIONS, sizeof(double));
    flagField = (unsigned int*) calloc(size,sizeof(int));
    Vels = (double*) calloc(size*3, sizeof(double));
    Temps = (double*) calloc(size*3, sizeof(double));
    
    /*initialization of fields*/
    initialiseFields(collideField_f, streamField_f, collideField_g, streamField_g, flagField, dim, velocityWall, filename);
    treatBoundary(collideField_f, collideField_g, flagField, Temps, velocityWall, dim, T_cold, T_warm);
    
    /*now start the calculation: */
    printf("=================================================================\n");
    printf("\nCOMPUTING...\n\n");
    for(int t = 1; t <= timesteps; t++){
        
        /*stream and put the results in collideField vector */
        double *swap_f=NULL;
        double *swap_g=NULL;
        doStreaming(collideField_f, streamField_f, collideField_g, streamField_g, flagField, dim);
        swap_f = collideField_f;
        swap_g = collideField_g;
        collideField_f = streamField_f;
        collideField_g = streamField_g;
        streamField_f = swap_f;
        streamField_g = swap_g;
        
        /*collide*/
        doCollision(collideField_f, collideField_g, flagField, &tau_f, &tau_g, dim, Vels, Temps, beta, T_cold, T_warm, gravity);

        /*apply boundary conditions*/
        treatBoundary(collideField_f, collideField_g, flagField, Temps, velocityWall, dim, T_cold, T_warm);

        /*write output and print progress to console*/
        if (t%timestepsPerPlotting==0){
            writeVtkOutput(Vels, Temps, flagField, argv[0], t, dim);
            printf("t = %4d/%d (%.1f%% completed)\n", t, timesteps, 100.0*t/timesteps);
        }
    }
    
    /*free up memory*/
    free(collideField_f);
    free(collideField_g);
    free(streamField_f);
    free(streamField_g);
    free(flagField);
    free(Vels);
    free(Temps);
    
    printf("\n=================================================================\n");
    printf("SIMULATION ENDED SUCCESFULLY");
	printf("\n=================================================================\n");
    return 0;
}

#endif

