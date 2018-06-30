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
    
    double *collideField=NULL;
    double *streamField=NULL;
    unsigned int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    unsigned int timesteps;
    unsigned int timestepsPerPlotting;
    double* vel;
    
    readParameters(&xlength,&tau,velocityWall,&timesteps,&timestepsPerPlotting,argc, argv);
    
    // Remove old results
    system("rm -rf Output");
    system("mkdir Output");
    
    /*memory allocation of required arrays*/
    int temp = (xlength+2) * (xlength+2) * (xlength+2);
    collideField = (double*) malloc(sizeof(double) * temp * NO_OF_LATTICE_DIMENSIONS);
    streamField = (double*) malloc(sizeof(double) * temp * NO_OF_LATTICE_DIMENSIONS);
    flagField = (unsigned int*) malloc(sizeof(int) * temp);
    vel = (double*)malloc(sizeof(double)*temp*3);
    
    /*initialization of Fields and printing*/
    initialiseFields(collideField,streamField,flagField,xlength, velocityWall);
  //  writeVtkOutput(collideField,flagField,argv[0],0,xlength, vel);
    
    for(int t = 1; t <= timesteps; t++){
        double *swap=NULL;
        doStreaming(collideField,streamField,flagField,xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        
        doCollision(collideField,flagField,&tau,xlength,vel);
        treatBoundary(collideField,flagField,velocityWall,xlength);
        
        if (t%timestepsPerPlotting==0){
            writeVtkOutput(collideField,flagField,argv[0],t,xlength, vel);
        }
    }
    /*free up memory*/
    free(collideField);
    free(streamField);
    free(flagField);
    free(vel);
    return 0;
}

#endif

