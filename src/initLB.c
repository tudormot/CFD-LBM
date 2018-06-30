#include "initLB.h"
#include "LBDefinitions.h"
#include "boundary.h"

int readParameters(int *xlength, double *tau, double *velocityWall, unsigned int *timesteps, unsigned int *timestepsPerPlotting, int argc, char *argv[]){

	/* check that there is only one command line parameter */
	if(argc !=2)
	{
		printf("ERROR! command line input not recognised. Was expecting filename..\n");
		return 0;
	}

	read_int( argv[1], "xlength", xlength);
	read_double( argv[1], "tau", tau);
	read_double( argv[1], "velocityWall", velocityWall);
	read_int( argv[1], "timesteps", (int *)timesteps);
	read_int( argv[1], "timestepsPerPlotting", (int *)timestepsPerPlotting);


  return 0;
}


void initialiseFields(double *collideField, double *streamField, unsigned int *flagField, int xlength ,const double * const wallVelocity){

    /* first initialise collideField and streamField: */
	/* also initialise the fluid cell flags in this for loop nest for efficiency*/
	for(int i = 1;i<=xlength;i++)
		for(int j = 1; j<=xlength;j++)
			for(int k = 1; k<=xlength;k++)
			{
				flagField[i*xlength*xlength + j*xlength + k] = FLUID;
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					collideField[NO_OF_LATTICE_DIMENSIONS*(i*xlength*xlength + j*xlength + k) + l]=LATTICEWEIGHTS[l];
					streamField[NO_OF_LATTICE_DIMENSIONS*(i*xlength*xlength + j*xlength + k) + l]=LATTICEWEIGHTS[l];
				}
			}

	/*Now initialise flag field on the boundaries*/

	for(int i = 0 ;i<=xlength+1;i++)
		for(int j = 0 ;j<=xlength+1;j++)
		{
			/*the lid of our cavity (top face) need to be set to moving wall*/
			flagField[i*xlength + j] = MOV_WALL;
			/*set bottom of cavity to no slip*/
			flagField[(xlength+1)*xlength*xlength+i*xlength + j] = NO_SLIP;
		}
	for(int k = 1;k<=xlength;k++)
		for(int i = 0;i<=xlength+1;i++)
		{
			flagField[k*xlength*xlength+ i] = NO_SLIP;
			flagField[k*xlength*xlength+(xlength+1)*xlength + i] = NO_SLIP;
		}
	for(int k = 1;k<=xlength;k++)
		for(int j = 1;j<=xlength;j++) //indices set so there is no overlapping, in case debugging is needed for loops could also overlap
		{
			flagField[k*xlength*xlength + j*xlength + 0] = NO_SLIP;
			flagField[k*xlength*xlength + j*xlength + (xlength +1)] = NO_SLIP;
		}

	/*finally set boundary values as well*/
	treatBoundary(collideField,flagField,wallVelocity,xlength);

}

