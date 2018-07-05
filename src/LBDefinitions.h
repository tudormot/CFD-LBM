#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

static const int LATTICEVELOCITIES[19][3] = {
    
    { 0, -1, -1},
    {-1,  0, -1},
    { 0,  0, -1},
    { 1,  0, -1},
    { 0,  1, -1},
    
    {-1, -1,  0},
    { 0, -1,  0},
    { 1, -1,  0},
    {-1,  0,  0},
    { 0,  0,  0},
    
    { 1,  0,  0},
    {-1,  1,  0},
    { 0,  1,  0},
    { 1,  1,  0},
    { 0, -1,  1},
    
    {-1,  0,  1},
    { 0,  0,  1},
    { 1,  0,  1},
    { 0,  1,  1}
    
};

static const double LATTICEWEIGHTS[19] = {1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0,
                                          1.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0, 12.0/36.0, 
                                          2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 
                                          1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0};
    
static const double C_S = 1.0/1.73205; // 1.73205 is estimated sqrt(3)

static const int R = 8.3144598; // Ideal gas constant
static const int D0 = 3;

#define NO_OF_DIMENSIONS (3)
#define NO_OF_LATTICE_DIMENSIONS (19)

/*definitions which have to do with flags*/
#define FLUID (0)
#define NO_SLIP (1)
#define MOV_WALL (2) //moving wall


typedef struct{
	int xlen;
	int ylen;
	int zlen;
}dimensions;

#endif
    
