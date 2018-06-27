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

#define c0  12.0/36.0 // TODO: Tudor will complain...
#define c1  2.0/36.0
#define c2  1.0/36.0

static const double LATTICEWEIGHTS[19] = {c2, c2, c1, c2, c2,
                                          c2, c1, c2, c1, c0, 
                                          c1, c2, c1, c2, c2, 
                                          c2, c1, c2, c2};
    
static const double C_S = 1.0/1.73205; // 1.73205 is estimated sqrt(3)
#endif
    
