#include "initLB.h"
#include "LBDefinitions.h"
#include "boundary.h"


static void specialInitFlags_NaturalConvection(unsigned int * flagField,dimensions dim);
static void specialInitFlags_Cavity(unsigned int * flagField,dimensions dim);



int readParameters(dimensions * dim, double *tau_f, double *tau_g, double *T_cold, double *T_warm, double *velocityWall, unsigned int *timesteps, unsigned int *timestepsPerPlotting, int argc, char *argv[], char *filename){


	/* check that there is only one command line parameter */
	if(argc !=2)
	{
		printf("ERROR! command line input not recognised. Was expecting filename..\n");
		return 0;
	}
	filename = argv[1];

	read_int( argv[1], "xlength", &(dim->xlen));
	read_int( argv[1], "ylength", &(dim->ylen));
	read_int( argv[1], "zlength", &(dim->zlen));
	read_double( argv[1], "tau_f", tau_f);
	read_double( argv[1], "tau_g", tau_f);
	read_double( argv[1], "T_cold", T_cold);
	read_double( argv[1], "T_warm", T_warm);
	read_double( argv[1], "xvelocityWall", &velocityWall[0]);
	read_double( argv[1], "yvelocityWall", &velocityWall[1]);
	read_double( argv[1], "zvelocityWall", &velocityWall[2]);
	read_int( argv[1], "timesteps", (int *)timesteps);
	read_int( argv[1], "timestepsPerPlotting", (int *)timestepsPerPlotting);


  return 0;
}


void initialiseFields(double *collideField_f, double *streamField_f,double *collideField_g, double *streamField_g,  unsigned int *flagField, dimensions dim ,const double * const wallVelocity, const char * const filename){

    /* first initialise collideField and streamField: */
	/* also initialise the fluid cell flags in this for loop nest for efficiency*/
	int xl = dim.xlen+2;           //variables used in the mapping between the spacial coordinates and the location in array
	int xlyl = xl*(dim.ylen+2);	   //see line above
	
	for(int i = 0;i<=dim.zlen+1;i++)
		for(int j = 0; j<=dim.ylen;j++)
			for(int k = 0; k<=dim.xlen+1;k++)
			{
				flagField[i*xlyl + j*xl + k] = IS_FLUID_BIT; //flag field  initialised to FLUID by default, reduces work to be done in specialsetflags functions
				for(int l = 0;l<NO_OF_LATTICE_DIMENSIONS;l++)
				{
					collideField_f[NO_OF_LATTICE_DIMENSIONS*(i*xlyl + j*xl + k) + l]=LATTICEWEIGHTS[l];
					streamField_f[NO_OF_LATTICE_DIMENSIONS*(i*xlyl + j*xl + k) + l]=LATTICEWEIGHTS[l];
					/*see document ../things_we_are_unsure.txt to why we are initialising g fields to 0: */
					{
						//actually there is no need to initialise g fields to 0 as they are initialised by calloc
					}
				}
			}

	/*Now initialise flag field on the boundaries. Based on input filename we are initialising flags differently:*/

	if(strcmp(filename,"cavity.dat")==0) //if name is cavity.dat
	{
		specialInitFlags_Cavity(flagField,dim);
	}
	else if (strcmp(filename,"convection.dat")==0)
	{
		specialInitFlags_NaturalConvection(flagField,dim);
	}
	else
	{
		printf("Error! no special_setflags function available for filename specified..\n");
	}




}


static void specialInitFlags_NaturalConvection(unsigned int * flagField,dimensions dim)
{
	int xl = dim.xlen+2;           //variables used in the mapping between the spacial coordinates and the location in array
	int xlyl = xl*(dim.ylen+2);	   //see line above

	/*setting opposite lateral faces which are adiabatic and free slip (to be able to compare with 2d navier stokes) */
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int i = 0;i<=dim.xlen+1;i++)
		{
			flagField[k*xlyl + i] = IS_FREESLIP_BIT | IS_NEUMANN_T_BIT;
			flagField[k*xlyl+(dim.ylen+1)*xl + i] = IS_FREESLIP_BIT | IS_NEUMANN_T_BIT;
		}

	/*setting opposite lateral faces of the cube with the heated walls and no slip velocity conditions*/
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int j = 0;j<=dim.ylen+1;j++) //indices are overlapping, however in this case since everything is non slip it does not matter
		{
			flagField[k*xlyl + j*xl + 0] = IS_NOSLIP_BIT | IS_DIRICHL_T_BIT | IS_COLD_WALL;
			flagField[k*xlyl + j*xl + (dim.xlen +1)] = IS_NOSLIP_BIT | IS_DIRICHL_T_BIT | IS_WARM_WALL;
		}
	/*setting flags for bottom and top faces of the cube*/
	for(int i = 0 ;i<=dim.ylen+1;i++)
		for(int j = 0 ;j<=dim.xlen+1;j++)
		{
			flagField[i*xl + j] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
			flagField[(dim.zlen+1)*xlyl+i*xl + j] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
		}
}

static void specialInitFlags_Cavity(unsigned int * flagField,dimensions dim)
{
	int xl = dim.xlen+2;           //variables used in the mapping between the spacial coordinates and the location in array
	int xlyl = xl*(dim.ylen+2);	   //see line above

	/*setting opposite lateral faces which are adiabatic and no slip*/
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int i = 0;i<=dim.xlen+1;i++)
		{
			flagField[k*xlyl + i] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
			flagField[k*xlyl+(dim.ylen+1)*xl + i] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
		}

	/*setting opposite lateral faces of the cube which are adiabatic and no slip*/
	for(int k = 0;k<=dim.zlen+1;k++)
		for(int j = 0;j<=dim.ylen+1;j++) //indices are overlapping, however in this case since everything is non slip it does not matter
		{
			flagField[k*xlyl + j*xl + 0] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
			flagField[k*xlyl + j*xl + (dim.xlen +1)] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
		}
	/*setting flags for bottom and top faces of the cube. Both are adiabatic, but the top face is a moving wall*/
	for(int i = 0 ;i<=dim.ylen+1;i++)
		for(int j = 0 ;j<=dim.xlen+1;j++)
		{
			flagField[i*xl + j] = IS_NOSLIP_BIT | IS_NEUMANN_T_BIT;
			flagField[(dim.zlen+1)*xlyl+i*xl + j] = IS_INFLOW_BIT | IS_NEUMANN_T_BIT;
		}
}
