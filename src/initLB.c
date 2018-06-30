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
	read_double( argv[1], "xvelocityWall", &velocityWall[0]);
	read_double( argv[1], "yvelocityWall", &velocityWall[1]);
	read_double( argv[1], "zvelocityWall", &velocityWall[2]);
	read_int( argv[1], "timesteps", (int *)timesteps);
	read_int( argv[1], "timestepsPerPlotting", (int *)timestepsPerPlotting);


  return 0;
}


void initialiseFields(double *collideField, double *streamField, unsigned int *flagField, int xlength ,const double * const wallVelocity){

    /* first initialise collideField and streamField: */
	/* also initialise the fluid cell flags in this for loop nest for efficiency*/
	int xl = xlength+2;
	int xl2 = xl*xl;
	
	for(int i = 0;i<=xlength+1;i++)
		for(int j = 0; j<=xlength+1;j++)
			for(int k = 0; k<=xlength+1;k++)
			{
				flagField[i*xl2 + j*xl + k] = FLUID;
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					collideField[NO_OF_LATTICE_DIMENSIONS*(i*xl2 + j*xl + k) + l]=LATTICEWEIGHTS[l];
					streamField[NO_OF_LATTICE_DIMENSIONS*(i*xl2 + j*xl + k) + l]=LATTICEWEIGHTS[l];
				}
			}

	/*Now initialise flag field on the boundaries*/

	
	for(int k = 0;k<=xlength+1;k++)
		for(int i = 0;i<=xlength+1;i++)
		{
			flagField[k*xl2 + i] = NO_SLIP;
			flagField[k*xl2+(xlength+1)*xl + i] = NO_SLIP;
		}
	for(int k = 0;k<=xlength+1;k++)
		for(int j = 0;j<=xlength+1;j++) //indices set so there is no overlapping, in case debugging is needed for loops could also overlap
		{
			flagField[k*xl2 + j*xl + 0] = NO_SLIP;
			flagField[k*xl2 + j*xl + (xlength +1)] = NO_SLIP;
		}
	for(int i = 0 ;i<=xlength+1;i++)
		for(int j = 0 ;j<=xlength+1;j++)
		{
			/*the lid of our cavity (top face) need to be set to moving wall*/
			flagField[i*xl + j] = NO_SLIP;
			/*set bottom of cavity to no slip*/
			flagField[(xlength+1)*xl2+i*xl + j] = MOV_WALL;
		}
	
}

