#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "LBDefinitions.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
dimensions * dim,                       /* reads domain size. Parameter name: "xlength" */
double *tau_f,                        /* relaxation parameter for f. Parameter name: "tau" */
double *tau_g,                        /* relaxation parameter for g. Parameter name: "tau" */
double *T_cold,						/* temperature of cold wall in convection case */
double *T_warm,						/* temperature of warm wall in convection case */
double *beta,						/* thermal expansion coefficient */
double *gravity,					/* gravity acceleration */
double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
unsigned int *timesteps,            /* number of timesteps. Parameter name: "timesteps" */
unsigned int *timestepsPerPlotting, /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
char *argv[],                        /* argv[1] shall contain the path to the config file */
char **filename                      /* initialised inside this function to argv[1]*/
);


/* initialises the particle distribution functions and the flagfield */
/* This function initialises the density fields as specified in worksheet 2.
 * It also initialises the flag field with the help of "special initialise functions,
 * which are chosen based on command line input"*/
void initialiseFields(double *collideField_f, double *streamField_f,double *collideField_g,
		double *streamField_g, unsigned int *flagField, dimensions dim,
		const double * const wallVelocity,const char * const filename);

#endif

