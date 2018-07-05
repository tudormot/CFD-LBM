#include "initLB.h"
#include "LBDefinitions.h"
#include "boundary.h"

int readParameters(dimensions * dim, double *tau_f, double *tau_g, double *velocityWall, unsigned int *timesteps, unsigned int *timestepsPerPlotting, int argc, char *argv[]){

	/* check that there is only one command line parameter */
	if(argc !=2)
	{
		printf("ERROR! command line input not recognised. Was expecting filename..\n");
		return 0;
	}

	read_int( argv[1], "xlength", &(dim->xlen));
	read_int( argv[1], "ylength", &(dim->ylen));
	read_int( argv[1], "zlength", &(dim->zlen));
	read_double( argv[1], "tau_f", tau_f);
	read_double( argv[1], "tau_g", tau_f);
	read_double( argv[1], "xvelocityWall", &velocityWall[0]);
	read_double( argv[1], "yvelocityWall", &velocityWall[1]);
	read_double( argv[1], "zvelocityWall", &velocityWall[2]);
	read_int( argv[1], "timesteps", (int *)timesteps);
	read_int( argv[1], "timestepsPerPlotting", (int *)timestepsPerPlotting);


  return 0;
}


void initialiseFields(double *collideField, double *streamField,  unsigned int *flagField, dimensions dim ,const double * const wallVelocity){

    /* first initialise collideField and streamField: */
	/* also initialise the fluid cell flags in this for loop nest for efficiency*/
	int xl = dim.xlen+2;           //variables used in the mapping between the spacial coordinates and the location in array
	int xlyl = xl*(dim.ylen+2);	   //see above
	
	for(int i = 0;i<=dim.zlen+1;i++)
		for(int j = 0; j<=dim.ylen;j++)
			for(int k = 0; k<=dim.xlen+1;k++)
			{
				flagField[i*xlyl + j*xl + k] = FLUID; //flag field  initialised to FLUID by default
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					collideField[NO_OF_LATTICE_DIMENSIONS*(i*xlyl + j*xl + k) + l]=LATTICEWEIGHTS[l];
					streamField[NO_OF_LATTICE_DIMENSIONS*(i*xlyl + j*xl + k) + l]=LATTICEWEIGHTS[l];
				}
			}

	/*Now initialise flag field on the boundaries*/
	initialiseFlags_NaturalConvection(flagField,dim);

	
#if 0 //old code
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int i = 0;i<=dim.xlen+1;i++)
		{
			flagField[k*xlyl + i] = NO_SLIP;
			flagField[k*xlyl+(dim.ylen+1)*xl + i] = NO_SLIP;
		}
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int j = 0;j<=dim.ylen+1;j++) //indices are overlapping, however in this case since everything is non slip it does not matter
		{
			flagField[k*xlyl + j*xl + 0] = NO_SLIP;
			flagField[k*xlyl + j*xl + (dim.xlen +1)] = NO_SLIP;
		}
	for(int i = 0 ;i<=dim.ylen+1;i++)
		for(int j = 0 ;j<=dim.xlen+1;j++)
		{
			/*the lid of our cavity (top face) need to be set to moving wall*/
			flagField[i*xl + j] = NO_SLIP;
			/*set bottom of cavity to no slip*/
			flagField[(dim.zlen+1)*xlyl+i*xl + j] = MOV_WALL;
		}
#endif
	
}

