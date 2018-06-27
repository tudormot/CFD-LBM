#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){
    
    double *collideField=NULL;
    double *streamField=NULL;
    int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    int timesteps;
    int timestepsPerPlotting;
    
    readParameters(
        &xlength,&tau,&velocityWall,&timesteps,&timestepsPerPlotting,argc, argv
    );
    
    /*memory allocation of required arrays and initialization:*/
    int temp = NO_OF_LATTICE_DIMENSIONS * (int)pow((double)(xlength+2),(double)NO_OF_DIMENSIONS); //(xlength + 2)^D
    collideField = (double*) malloc(sizeof(double) * temp * NO_OF_LATTICE_DIMENSIONS);
    streamField = (double*) malloc(sizeof(double) * temp * NO_OF_LATTICE_DIMENSIONS);
    flagField = (int*) malloc(sizeof(int) * temp);
    initialiseFields(collideField,streamField,flagField,xlength);
    
    for(int t = 0; t < timesteps; t++){
        double *swap=NULL;
        doStreaming(collideField,streamField,flagField,xlength);
        swap = collideField;
        collideField = streamField;
        streamField = swap;
        
        doCollision(collideField,flagField,&tau,xlength);
        treatBoundary(collideField,flagField,velocityWall,xlength);
        
        if (t%timestepsPerPlotting==0){
            writeVtkOutput(collideField,flagField,argv,t,xlength);
        }
    }
    /*free up memory*/
    free(collideField);
    free(streamField);
    free(flagField);
    
    return 0;
}

#endif

