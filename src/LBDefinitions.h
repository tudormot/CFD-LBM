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
    
}

double c0 = 12.0/36.0;
double c1 =  2.0/36.0;
double c2 =  1.0/36.0;

static const double LATTICEWEIGHTS[19] = {c2, c2, c1, c2, c2, 
                                          c2, c1, c2, c1, c0, 
                                          c1, c2, c1, c2, c2, 
                                          c2, c1, c2, c2};
    
static const double C_S = 1.0/sqrt(3.0);
#endif
    
